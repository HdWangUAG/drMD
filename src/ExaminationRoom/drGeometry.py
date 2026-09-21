## BASIC PYTHON LIBRARIES
import os
from os import path as p
import argparse
import json

## NUMERICAL LIBRARIES
import numpy as np
import pandas as pd
import yaml

## CLEAN CODE
from typing import Dict, List, Optional, Sequence, Tuple

########################################################################################################
"""
Geometric observables of a drMD trajectory: distances, centre-of-mass distances, minimum heavy-atom
distances and angles between groups of atoms, measured with the topology that was actually simulated.

This is the counterpart of src/ExaminationRoom/drPLIP.py. PLIP tells you which interactions are
present; this tells you how far apart a chosen pair of atoms or groups is, frame by frame, using the
hydrogens and tautomers of the simulation rather than a re-protonated structure. The same
measurements define metadynamics collective variables, so a measurement file doubles as the
justification for a CV range.

A measurement file is YAML (or JSON):

    chainMap: "A:1-299:83,B:300-598:83,C:599-675:0,D:676-752:0"
    measurements:
      - name: C1-His285_NE2
        type: distance
        a: {chain: A, resId: 285, atom: NE2}
        b: {chain: C, resId: 36, atom: C1}
      - name: tail-pocketBottom
        type: comDistance
        a: {chain: C, resId: 36, atoms: [C10, C11, C12]}
        b: {chain: A, resIds: [137, 140, 141, 146, 189, 199], atoms: sidechain}
      - name: Met197-tail_closest
        type: minDistance
        a: {chain: A, resId: 197, atoms: sidechain}
        b: {chain: C, resId: 36, atoms: [C10, C11, C12]}
      - name: His285_chi1
        type: torsion
        a: {chain: A, resId: 285, atom: N}
        b: {chain: A, resId: 285, atom: CA}
        c: {chain: A, resId: 285, atom: CB}
        d: {chain: A, resId: 285, atom: CG}
      - name: thioester-His285
        type: angle
        a: {chain: C, resId: 36, atom: O1}
        b: {chain: C, resId: 36, atom: C1}
        c: {chain: A, resId: 285, atom: NE2}
      - name: C1-nearestWater
        type: nearestSolvent
        a: {chain: C, resId: 36, atom: C1}
      - name: waterBridge_C1-His285
        type: solventBridge
        a: {chain: C, resId: 36, atom: C1}
        b: {chain: A, resId: 285, atom: NE2}
        cutoffA: 4.0
        cutoffB: 3.5

Selections: `chain` plus `resId`/`resIds`, and either `atom`/`atoms` (names), or
`atoms: sidechain | heavy | backbone | all`. Residue numbers follow `chainMap` when one is given,
otherwise the numbering of the topology PDB.

`minDistance` is the closest approach of any heavy atom of `a` to any heavy atom of `b`, per frame - the
number a structural reviewer means by "the residue is within 4 Å of the ligand". A centroid distance
(`comDistance`) can read 8 Å while the side chain is in van der Waals contact, so it is not a substitute.
Hydrogens are excluded by element whatever the selection says: `atoms: all` and `atoms: heavy` give the
same answer, and a hydrogen named explicitly is refused rather than dropped without a word. Every pair is
evaluated with the minimum-image convention, so a contact between two chains that straddle a periodic
boundary is measured correctly. The number of heavy-atom pairs is capped (`maxPairs` on the measurement,
default MIN_DISTANCE_MAX_PAIRS) so that two whole chains are refused rather than ground through.

The solvent types (nearestSolvent, nearestSolventTo, solventBridge, bridgeAngle) look at the water
oxygens of the simulated box, with the minimum-image convention, because the water that matters is
often the one that has just diffused across a periodic boundary. They count *populations* of a
geometry - how many waters sit in a bridging position, how often a water is in reach - and are not
rates and not barriers.

Outputs (in --outDir): geometry_series.csv (one row per frame), geometry_summary.csv / .md
(mean, sd, min, max and the fraction of frames below each --cutoff) and geometry_run.json.
"""
########################################################################################################

BACKBONE = {"N", "CA", "C", "O", "OXT"}

## residue names counted as solvent on top of mdtraj's own residue.is_water, which does not know the
## names that tleap and CHARMM-GUI write for a three-site water
SOLVENT_RESIDUES = {"HOH", "WAT", "T3P", "TIP3", "SOL"}
## a water is represented by its oxygen: hydrogens wag, the oxygen is where the lone pairs are
SOLVENT_OXYGEN_NAMES = {"O", "OW", "OH2"}
SOLVENT_TYPES = ("nearestSolvent", "nearestSolventTo", "solventBridge", "bridgeAngle")

## minDistance evaluates every heavy-atom pair between its two groups. A residue against a ligand is a few
## hundred pairs and a ligand against a whole chain ~100,000; two whole chains run to millions, which is
## a contact map rather than a measurement, so the type refuses above this many pairs unless the
## measurement raises `maxPairs` deliberately. Pairs are evaluated in blocks so that the (frames, pairs)
## matrix of a chunk stays small however large the groups are.
MIN_DISTANCE_MAX_PAIRS = 2_000_000
MIN_DISTANCE_PAIR_BLOCK = 50_000


def parse_chain_map(chainMap: Optional[str]) -> List[Tuple[str, int, int, int]]:
    """'A:1-299:83,C:599-675:0' -> [(chainId, firstResidueIndex1, lastResidueIndex1, numberingOffset)]"""
    if not chainMap:
        return []
    entries = []
    for item in chainMap.split(","):
        chainId, span, offset = item.strip().split(":")
        first, last = span.split("-")
        entries.append((chainId, int(first), int(last), int(offset)))
    return entries


def is_bulk_residue(residue, solventResidues: Optional[Sequence[str]] = None) -> bool:
    """Water, or a monatomic ion: the bath, not something a measurement selects by residue number."""
    names = SOLVENT_RESIDUES if solventResidues is None else {n.upper() for n in solventResidues}
    if getattr(residue, "is_water", False) or residue.name.upper() in names:
        return True
    return residue.n_atoms == 1


def residue_labels(topology, chainMap: Optional[str],
                   solventResidues: Optional[Sequence[str]] = None) -> Dict[Tuple[str, int], int]:
    """(chainId, residueNumber) -> residue index in the topology.

    Residues covered by the chainMap are labelled first and must be unique - a clash there means the
    chainMap does not separate the chains, which is a mistake worth stopping for. Everything else is
    labelled by its PDB numbering only when that is unambiguous: a solvated box holds tens of thousands
    of waters whose resSeq wraps round at 9999, so insisting on uniqueness there would refuse every
    explicit-solvent trajectory. Those residues simply cannot be named in a measurement (the solvent
    types find them by element instead).
    """
    mapping = parse_chain_map(chainMap)
    labels: Dict[Tuple[str, int], int] = {}
    unmapped: List[Tuple[Tuple[str, int], object]] = []
    for residue in topology.residues:
        label = None
        for chainId, first, last, offset in mapping:
            if first <= residue.index + 1 <= last:
                label = (chainId, residue.index + 1 - first + 1 + offset)
                break
        if label is None:
            chain = residue.chain
            chainId = getattr(chain, "chain_id", None) or chr(ord("A") + chain.index % 26)
            unmapped.append(((chainId, residue.resSeq), residue))
            continue
        if label in labels:
            raise ValueError(f"residue {label[0]}{label[1]} is not unique; give a --chainMap that separates the chains")
        labels[label] = residue.index

    clashing = {label for label, _ in unmapped if label in labels}
    seen = set()
    for label, _ in unmapped:
        if label in seen:
            clashing.add(label)
        seen.add(label)
    for label, residue in unmapped:
        if label not in clashing:
            labels[label] = residue.index
        elif not is_bulk_residue(residue, solventResidues):
            raise ValueError(f"residue {label[0]}{label[1]} is not unique; give a --chainMap that separates the chains")
    return labels


def select_atoms(topology, labels: Dict[Tuple[str, int], int], selection: dict) -> List[int]:
    """Resolve one selection dictionary to a list of atom indices."""
    chainId = selection.get("chain")
    resIds = selection.get("resIds", [selection["resId"]] if "resId" in selection else None)
    if chainId is None or resIds is None:
        raise ValueError(f"selection needs 'chain' and 'resId'/'resIds': {selection}")
    wanted = selection.get("atoms", selection.get("atom", "heavy"))
    names = None
    if isinstance(wanted, str) and wanted not in ("sidechain", "heavy", "backbone", "all"):
        names = {wanted}
    elif isinstance(wanted, (list, tuple)):
        names = set(wanted)

    indices = []
    for resId in resIds:
        key = (chainId, int(resId))
        if key not in labels:
            raise KeyError(f"residue {chainId}{resId} is not in the topology")
        ## topology.residue(i) is a direct lookup; list(topology.residues) would walk 50,000 waters
        ## on every selection of every chunk of a solvated trajectory
        residue = topology.residue(labels[key])
        for atom in residue.atoms:
            isHydrogen = atom.element is not None and atom.element.symbol == "H"
            if names is not None:
                if atom.name in names:
                    indices.append(atom.index)
            elif wanted == "all":
                indices.append(atom.index)
            elif wanted == "heavy" and not isHydrogen:
                indices.append(atom.index)
            elif wanted == "backbone" and atom.name in BACKBONE:
                indices.append(atom.index)
            elif wanted == "sidechain" and not isHydrogen and atom.name not in BACKBONE:
                indices.append(atom.index)
    if not indices:
        raise ValueError(f"selection matched no atoms: {selection}")
    if names is not None:
        missing = names - {topology.atom(i).name for i in indices}
        if missing:
            raise ValueError(f"atoms {sorted(missing)} not found for selection {selection}")
    return indices


########################################################################################################
def centroid(xyz: np.ndarray, masses: Optional[np.ndarray]) -> np.ndarray:
    """xyz: (frames, atoms, 3). Mass-weighted if masses are given, else geometric."""
    if masses is None:
        return xyz.mean(axis=1)
    weights = masses / masses.sum()
    return (xyz * weights[None, :, None]).sum(axis=1)


def solvent_oxygen_indices(topology, solventResidues: Optional[Sequence[str]] = None) -> np.ndarray:
    """Indices of the water oxygens: mdtraj's is_water, widened by a set of residue names."""
    names = SOLVENT_RESIDUES if solventResidues is None else {n.upper() for n in solventResidues}
    indices = []
    for residue in topology.residues:
        if not (getattr(residue, "is_water", False) or residue.name.upper() in names):
            continue
        for atom in residue.atoms:
            isOxygen = atom.element is not None and atom.element.symbol == "O"
            if atom.name.upper() in SOLVENT_OXYGEN_NAMES and isOxygen:
                indices.append(atom.index)
    return np.array(sorted(indices), dtype=int)


def probe_trajectory(traj, atomIndices: np.ndarray, points: Sequence[np.ndarray]):
    """A small trajectory of the given atoms plus one pseudo-atom per point, carrying the same box.

    Centroids are not atoms, so mdtraj cannot measure to them directly; and copying the whole 160,000-atom
    frame to add three sites would cost more than the measurement. Building this stand-in keeps every
    solvent distance inside mdtraj's minimum-image code, which is the only thing here that gets a triclinic
    box right. Coordinates are nanometres, as everywhere in mdtraj.
    """
    import mdtraj as md

    atomIndices = np.asarray(atomIndices, dtype=int)
    stacked = [traj.xyz[:, atomIndices, :]] + [np.asarray(point, dtype=np.float32).reshape(traj.n_frames, 1, 3)
                                               for point in points]
    xyz = np.concatenate(stacked, axis=1)
    probe = md.Trajectory(xyz=xyz, topology=None,
                          unitcell_lengths=traj.unitcell_lengths, unitcell_angles=traj.unitcell_angles)
    atomSlots = np.arange(len(atomIndices))
    pointSlots = [len(atomIndices) + i for i in range(len(points))]
    return probe, atomSlots, pointSlots


def per_frame_distances(traj, pairs: np.ndarray, periodic: bool = True, block: int = 256) -> np.ndarray:
    """Distance (Angstrom) of one pair per frame, where the pair may differ from frame to frame.

    md.compute_distances evaluates EVERY pair in EVERY frame and returns (nFrames, nPairs). Handing it
    one pair per frame therefore gives an (nFrames, nFrames) matrix whose wanted values are the diagonal;
    taking [:, 0] would silently report frame 0's water in every frame, which is exactly the shape of bug
    that survives a review because the numbers look plausible. Blocking keeps that square matrix small.
    """
    import mdtraj as md

    pairs = np.asarray(pairs, dtype=int)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError(f"per_frame_distances needs pairs of shape (nFrames, 2), got {pairs.shape}")
    if pairs.shape[0] != traj.n_frames:
        raise ValueError(f"per_frame_distances needs one pair per frame: {pairs.shape[0]} pairs, {traj.n_frames} frames")
    out = np.empty(traj.n_frames, dtype=float)
    for start in range(0, traj.n_frames, block):
        stop = min(start + block, traj.n_frames)
        square = md.compute_distances(traj[start:stop], pairs[start:stop], periodic=periodic)
        out[start:stop] = np.diagonal(square)
    return out * 10.0


def per_frame_angles(traj, triplets: np.ndarray, periodic: bool = True, block: int = 256) -> np.ndarray:
    """Angle (degrees) of one triplet per frame, vertex in the middle. Diagonal, for the reason above."""
    import mdtraj as md

    triplets = np.asarray(triplets, dtype=int)
    if triplets.ndim != 2 or triplets.shape[1] != 3:
        raise ValueError(f"per_frame_angles needs triplets of shape (nFrames, 3), got {triplets.shape}")
    if triplets.shape[0] != traj.n_frames:
        raise ValueError(f"per_frame_angles needs one triplet per frame: {triplets.shape[0]} triplets, "
                         f"{traj.n_frames} frames")
    out = np.empty(traj.n_frames, dtype=float)
    for start in range(0, traj.n_frames, block):
        stop = min(start + block, traj.n_frames)
        square = md.compute_angles(traj[start:stop], triplets[start:stop], periodic=periodic)
        out[start:stop] = np.degrees(np.diagonal(square))
    return out


def heavy_atom_indices(topology, indices: Sequence[int], selection: dict, name: str) -> np.ndarray:
    """The non-hydrogen atoms of a resolved selection, for the heavy-atom contact types.

    Hydrogens are recognised by element rather than by name, because ligand hydrogens do not follow
    the protein naming conventions. A class selection (`atoms: all`) simply loses its hydrogens - the
    measurement is defined as a heavy-atom distance, so `all` and `heavy` mean the same thing here. An
    atom the user named explicitly is different: quietly ignoring `atom: HB1` would return a number that
    is not the one asked for, and that is the kind of mistake that survives a review, so it is refused.
    """
    wanted = selection.get("atoms", selection.get("atom", "heavy"))
    namedExplicitly = isinstance(wanted, (list, tuple)) or (
        isinstance(wanted, str) and wanted not in ("sidechain", "heavy", "backbone", "all"))
    heavy, hydrogens = [], []
    for index in indices:
        atom = topology.atom(index)
        if atom.element is not None and atom.element.symbol == "H":
            hydrogens.append(atom.name)
        else:
            heavy.append(index)
    if hydrogens and namedExplicitly:
        raise ValueError(f"{name}: minDistance is a heavy-atom distance, but the selection names the "
                         f"hydrogen(s) {sorted(set(hydrogens))}: {selection}")
    if not heavy:
        raise ValueError(f"{name}: selection has no heavy atoms: {selection}")
    return np.asarray(heavy, dtype=int)


def measure_min_distance(traj, indicesA: np.ndarray, indicesB: np.ndarray, name: str,
                         maxPairs: int = MIN_DISTANCE_MAX_PAIRS) -> np.ndarray:
    """Closest approach (Angstrom) between any atom of A and any atom of B, per frame.

    Every pair goes through mdtraj's minimum-image code, like the solvent types: a TE chain and an ACP
    chain can sit on opposite sides of a periodic boundary, and a naive norm would report the distance
    to the wrong image. An exact all-pairs minimum is used rather than a cell list because, at the sizes
    the cap allows, mdtraj evaluates a million pairs over a hundred frames in well under a second; the
    cap exists so that an accidental chain-against-chain selection fails at once instead of quietly
    turning a two-minute analysis into an afternoon.
    """
    import mdtraj as md

    indicesA = np.asarray(indicesA, dtype=int)
    indicesB = np.asarray(indicesB, dtype=int)
    nPairs = indicesA.size * indicesB.size
    if nPairs > maxPairs:
        raise ValueError(f"{name}: minDistance between {indicesA.size} and {indicesB.size} heavy atoms is "
                         f"{nPairs} pairs, above the limit of {maxPairs}; narrow the selections (a pocket "
                         f"rather than a chain) or set 'maxPairs' on the measurement to insist")
    pairs = np.array(np.meshgrid(indicesA, indicesB, indexing="ij")).reshape(2, -1).T
    closest = np.full(traj.n_frames, np.inf)
    for start in range(0, nPairs, MIN_DISTANCE_PAIR_BLOCK):
        block = md.compute_distances(traj, pairs[start:start + MIN_DISTANCE_PAIR_BLOCK], periodic=True)
        np.minimum(closest, block.min(axis=1), out=closest)
    return closest * 10.0


def measure_solvent(traj, topology, groups: Sequence[np.ndarray], measurement: dict,
                    solventIndices: Optional[np.ndarray] = None,
                    solventResidues: Optional[Sequence[str]] = None) -> np.ndarray:
    """The solvent measurement types. `groups` are the centroids of a/b/c in Angstrom.

    Every distance here is minimum-image (periodic=True): the bridging water is frequently an image of a
    water written on the far side of the box, and a non-periodic distance would quietly miss it.
    """
    import mdtraj as md

    kind = measurement["type"]
    name = measurement.get("name", kind)
    oxygens = solvent_oxygen_indices(topology, solventResidues) if solventIndices is None else np.asarray(solventIndices)
    if oxygens.size == 0:
        raise ValueError(f"{name}: no solvent oxygens in the topology; is this an implicit-solvent run? "
                         f"(--solventResidues sets the residue names that count as solvent)")
    needed = {"nearestSolvent": 1, "nearestSolventTo": 2, "solventBridge": 2, "bridgeAngle": 3}[kind]
    if len(groups) != needed:
        raise ValueError(f"{name}: '{kind}' needs {needed} selection(s) "
                         f"({', '.join(['a', 'b', 'c'][:needed])}), got {len(groups)}")

    probe, oxygenSlots, pointSlots = probe_trajectory(traj, oxygens, [group / 10.0 for group in groups])
    ## distance from selection 'a' to every water oxygen, one column per water
    pairsA = np.column_stack([np.full(oxygenSlots.size, pointSlots[0]), oxygenSlots])
    toA = md.compute_distances(probe, pairsA, periodic=True) * 10.0

    if kind == "nearestSolvent":
        return toA.min(axis=1)
    if kind == "solventBridge":
        cutoffA = float(measurement.get("cutoffA", 4.0))
        cutoffB = float(measurement.get("cutoffB", 3.5))
        pairsB = np.column_stack([np.full(oxygenSlots.size, pointSlots[1]), oxygenSlots])
        toB = md.compute_distances(probe, pairsB, periodic=True) * 10.0
        return ((toA < cutoffA) & (toB < cutoffB)).sum(axis=1).astype(float)

    ## the remaining types follow one water, and which water that is changes from frame to frame
    nearest = oxygenSlots[toA.argmin(axis=1)]
    if kind == "nearestSolventTo":
        pairs = np.column_stack([nearest, np.full(traj.n_frames, pointSlots[1])])
        return per_frame_distances(probe, pairs)
    if kind == "bridgeAngle":
        triplets = np.column_stack([nearest, np.full(traj.n_frames, pointSlots[1]),
                                    np.full(traj.n_frames, pointSlots[2])])
        return per_frame_angles(probe, triplets)
    raise ValueError(f"unknown solvent measurement type '{kind}'")


def point_probe(traj, points: Sequence[np.ndarray]):
    """A probe trajectory holding just these per-frame points, so mdtraj's minimum-image code applies.

    A DCD holds wrapped coordinates. When a partner chain drifts across a box face it comes back on
    the other side, and a plain subtraction then reports a separation of order the box vector - tens
    of Angstrom - for a complex that never came apart. Routing every distance, angle and torsion
    through mdtraj keeps them on the minimum image, which is also the only code here that handles a
    triclinic (octahedral) box correctly. Points arrive in Angstrom and are converted, because this
    module works in Angstrom while mdtraj works in nanometres.
    """
    probe, _, pointSlots = probe_trajectory(traj, np.array([], dtype=int),
                                            [point / 10.0 for point in points])
    return probe, pointSlots


def point_distances(traj, points: Sequence[np.ndarray]) -> np.ndarray:
    """Minimum-image distance in Angstrom between two per-frame points."""
    import mdtraj as md

    probe, slots = point_probe(traj, points)
    return md.compute_distances(probe, np.array([slots[:2]]), periodic=True)[:, 0] * 10.0


def point_angles(traj, points: Sequence[np.ndarray]) -> np.ndarray:
    """Minimum-image angle in degrees; the middle point is the vertex."""
    import mdtraj as md

    probe, slots = point_probe(traj, points)
    return np.degrees(md.compute_angles(probe, np.array([slots[:3]]), periodic=True)[:, 0])


def point_torsions(traj, points: Sequence[np.ndarray]) -> np.ndarray:
    """Minimum-image torsion in degrees, in (-180, 180].

    mdtraj's dihedral follows the same sign convention as the CustomTorsionForce restraint in
    drRestraints, so a chi1 held by a restraint can be compared directly with the value measured
    after the restraint is released.
    """
    import mdtraj as md

    probe, slots = point_probe(traj, points)
    return np.degrees(md.compute_dihedrals(probe, np.array([slots[:4]]), periodic=True)[:, 0])


def measure(traj, topology, labels, measurement: dict, massWeighted: bool = True,
            solventIndices: Optional[np.ndarray] = None,
            solventResidues: Optional[Sequence[str]] = None) -> np.ndarray:
    """Return the measurement for every frame, in Angstrom (distances), degrees (angles) or a count."""
    kind = measurement.get("type", "distance")
    groups = []
    for key in ("a", "b", "c", "d"):
        if key in measurement:
            indices = select_atoms(topology, labels, measurement[key])
            masses = np.array([topology.atom(i).element.mass for i in indices]) if massWeighted else None
            groups.append(centroid(traj.xyz[:, indices, :] * 10.0, masses))
    if kind in SOLVENT_TYPES:
        return measure_solvent(traj, topology, groups, measurement, solventIndices, solventResidues)
    if kind == "minDistance":
        if len(groups) != 2:
            raise ValueError(f"{measurement.get('name')}: a minDistance needs selections 'a' and 'b'")
        name = measurement.get("name", kind)
        heavy = [heavy_atom_indices(topology, select_atoms(topology, labels, measurement[key]), measurement[key], name)
                 for key in ("a", "b")]
        return measure_min_distance(traj, heavy[0], heavy[1], name,
                                    maxPairs=int(measurement.get("maxPairs", MIN_DISTANCE_MAX_PAIRS)))
    if kind in ("distance", "comDistance"):
        if len(groups) != 2:
            raise ValueError(f"{measurement.get('name')}: a distance needs selections 'a' and 'b'")
        if kind == "distance":
            for key in ("a", "b"):
                if len(select_atoms(topology, labels, measurement[key])) != 1:
                    raise ValueError(f"{measurement.get('name')}: 'distance' needs one atom per selection; "
                                     f"use type 'comDistance' (centroids) or 'minDistance' (closest "
                                     f"heavy atoms) for a group")
        return point_distances(traj, groups)
    if kind == "angle":
        if len(groups) != 3:
            raise ValueError(f"{measurement.get('name')}: an angle needs selections 'a', 'b' and 'c' (b is the vertex)")
        return point_angles(traj, groups)
    if kind == "torsion":
        if len(groups) != 4:
            raise ValueError(f"{measurement.get('name')}: a torsion needs selections 'a', 'b', 'c' and 'd', "
                             f"in the order the dihedral runs")
        return point_torsions(traj, groups)
    raise ValueError(f"unknown measurement type '{kind}'")


def circular_statistics(degrees: np.ndarray) -> Tuple[float, float]:
    """Mean and spread of an angular series, so that -179 and +179 are two degrees apart.

    The mean is the direction of the resultant vector, and the spread is the circular standard
    deviation sqrt(-2 ln R), which tends to the ordinary standard deviation for a tight cluster.
    """
    radians = np.radians(degrees)
    resultant = np.hypot(np.sin(radians).mean(), np.cos(radians).mean())
    mean = np.degrees(np.arctan2(np.sin(radians).mean(), np.cos(radians).mean()))
    spread = np.degrees(np.sqrt(-2.0 * np.log(resultant))) if resultant > 0 else float("nan")
    return float(mean), float(spread)


def summarise(series: pd.DataFrame, cutoffs: Sequence[float],
              circularNames: Sequence[str] = ()) -> pd.DataFrame:
    rows = []
    for name in series.columns:
        if name in ("frame", "timeNs"):
            continue
        values = series[name].to_numpy()
        if name in circularNames:
            ## torsions wrap, so a linear mean of a series straddling 180 degrees is meaningless.
            ## min, max and the percentiles are reported about the circular mean for the same reason.
            mean, spread = circular_statistics(values)
            centred = (values - mean + 180) % 360 - 180
            row = {"name": name, "mean": mean, "sd": spread,
                   "min": mean + centred.min(), "max": mean + centred.max(),
                   "p5": mean + np.percentile(centred, 5), "p95": mean + np.percentile(centred, 95)}
        else:
            row = {"name": name, "mean": values.mean(), "sd": values.std(), "min": values.min(),
                   "max": values.max(), "p5": np.percentile(values, 5), "p95": np.percentile(values, 95)}
        for cutoff in cutoffs:
            row[f"frac_lt_{cutoff:g}"] = float("nan") if name in circularNames else float((values < cutoff).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def summary_markdown(summary: pd.DataFrame, nFrames: int, cutoffs: Sequence[float]) -> str:
    cutoffColumns = [f"frac_lt_{c:g}" for c in cutoffs]
    header = "| measurement | mean | sd | min | max | " + " | ".join(cutoffColumns) + " |\n"
    header += "|---|---|---|---|---|" + "---|" * len(cutoffColumns) + "\n"
    lines = [f"Geometry over {nFrames} frames (distances in Å, angles and torsions in degrees, "
             f"solventBridge in waters; torsions use circular statistics).\n\n", header]
    for _, r in summary.iterrows():
        cells = " | ".join("-" if pd.isna(r[c]) else f"{r[c]:.2f}" for c in cutoffColumns)
        lines.append(f"| {r['name']} | {r['mean']:.2f} | {r['sd']:.2f} | {r['min']:.2f} | {r['max']:.2f} | {cells} |\n")
    return "".join(lines)


########################################################################################################
def main() -> None:
    parser = argparse.ArgumentParser(description="Measure distances, COM distances, minimum heavy-atom distances, angles, torsions "
                                                 "and water bridges over a drMD trajectory.")
    parser.add_argument("--pdb", required=True, help="topology PDB")
    parser.add_argument("--trajectory", default=None, help="DCD trajectory; omit to measure the PDB alone")
    parser.add_argument("--measurements", required=True, help="YAML/JSON file describing the measurements")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--frameTimeNs", type=float, default=None,
                        help="time between the USED frames in ns, i.e. logInterval x stride. Give it "
                             "explicitly: a DCD may carry no per-frame time, and then there is none to infer")
    parser.add_argument("--chainMap", default=None, help="overrides the chainMap of the measurement file")
    parser.add_argument("--cutoffs", nargs="*", type=float, default=[3.5, 4.0],
                        help="report the fraction of frames below each of these values")
    parser.add_argument("--geometricCentre", action="store_true", help="use geometric instead of mass-weighted centroids")
    parser.add_argument("--chunk", type=int, default=100,
                        help="frames held in memory at once; a solvated trajectory does not fit whole")
    parser.add_argument("--solventResidues", nargs="*", default=sorted(SOLVENT_RESIDUES),
                        help="residue names counted as solvent, on top of mdtraj's is_water")
    parser.add_argument("--outDir", default="geometry")
    args = parser.parse_args()

    import mdtraj as md

    with open(args.measurements) as fh:
        spec = yaml.safe_load(fh)
    measurements = spec["measurements"] if isinstance(spec, dict) else spec
    chainMap = args.chainMap or (spec.get("chainMap") if isinstance(spec, dict) else None)

    topology = md.load(args.pdb).topology
    labels = residue_labels(topology, chainMap, args.solventResidues)
    ## scanning 50,000 residues once, rather than once per measurement per chunk
    solventIndices = (solvent_oxygen_indices(topology, args.solventResidues)
                      if any(m.get("type") in SOLVENT_TYPES for m in measurements) else None)
    if solventIndices is not None:
        print(f"solvent: {solventIndices.size} water oxygens", flush=True)

    columns = {measurement["name"]: [] for measurement in measurements}
    times: List[np.ndarray] = []
    if args.trajectory:
        ## streamed, not md.load: a production box of ~160,000 atoms is several GB on disk and the
        ## whole-trajectory load would be the only thing in this module that cannot finish
        for chunk in md.iterload(args.trajectory, top=topology, chunk=max(1, args.chunk), stride=args.stride):
            for measurement in measurements:
                columns[measurement["name"]].append(
                    measure(chunk, topology, labels, measurement, massWeighted=not args.geometricCentre,
                            solventIndices=solventIndices, solventResidues=args.solventResidues))
            times.append(np.asarray(chunk.time, dtype=float))
            print(f"  ...{sum(len(t) for t in times)} frames", flush=True)
    else:
        frame = md.load(args.pdb)
        for measurement in measurements:
            columns[measurement["name"]].append(
                measure(frame, topology, labels, measurement, massWeighted=not args.geometricCentre,
                        solventIndices=solventIndices, solventResidues=args.solventResidues))
        times.append(np.asarray(frame.time, dtype=float))

    frameTimes = np.concatenate(times) if times else np.zeros(0)
    nFrames = int(frameTimes.size)
    ## The time axis has to be right or event lifetimes are read off a wrong clock. Two traps:
    ## OpenMM's DCDs here carry no usable time (every stamp identical, so consecutive differences
    ## are zero), and --stride means the analysed series is coarser than the file by that factor.
    ## An unknown spacing is recorded as unknown rather than silently becoming zero.
    dt, dtSource = args.frameTimeNs, "given"
    if dt is None:
        spacing = np.diff(frameTimes) if frameTimes.size > 1 else np.zeros(0)
        usable = spacing[spacing > 0]
        if usable.size:
            dt = float(np.median(usable)) / 1000.0 * args.stride
            dtSource = f"trajectory median spacing x stride {args.stride}"
        else:
            dt, dtSource = None, "unknown: the trajectory carries no per-frame time"
            print("  no usable per-frame time in the trajectory; timeNs left as NaN. "
                  "Pass --frameTimeNs (logInterval x stride) if you need a time axis", flush=True)
    timeColumn = np.arange(nFrames) * dt if dt is not None else np.full(nFrames, np.nan)
    series = pd.DataFrame({"frame": np.arange(nFrames), "timeNs": timeColumn})
    for measurement in measurements:
        name = measurement["name"]
        series[name] = np.concatenate(columns[name])
        print(f"  {name}: mean {series[name].mean():.2f}, range [{series[name].min():.2f}, {series[name].max():.2f}]", flush=True)

    os.makedirs(args.outDir, exist_ok=True)
    series.to_csv(p.join(args.outDir, "geometry_series.csv"), index=False)
    torsionNames = [m["name"] for m in measurements if m.get("type") == "torsion"]
    summary = summarise(series, args.cutoffs, circularNames=torsionNames)
    summary.to_csv(p.join(args.outDir, "geometry_summary.csv"), index=False)
    with open(p.join(args.outDir, "geometry_summary.md"), "w") as fh:
        fh.write(summary_markdown(summary, nFrames, args.cutoffs))
    with open(p.join(args.outDir, "geometry_run.json"), "w") as fh:
        json.dump({"pdb": args.pdb, "trajectory": args.trajectory, "stride": args.stride,
                   "chunk": args.chunk, "nFrames": nFrames, "frameTimeNs": dt,
                   "frameTimeSource": dtSource, "chainMap": chainMap,
                   "cutoffs": list(args.cutoffs), "massWeighted": not args.geometricCentre,
                   "solventResidues": list(args.solventResidues),
                   "nSolventOxygens": int(solventIndices.size) if solventIndices is not None else 0,
                   "measurements": measurements}, fh, indent=2)
    print(f"done: {len(measurements)} measurements over {nFrames} frames -> {args.outDir}")


if __name__ == "__main__":
    main()
