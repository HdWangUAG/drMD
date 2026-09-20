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
Geometric observables of a drMD trajectory: distances, centre-of-mass distances and angles between
groups of atoms, measured with the topology that was actually simulated.

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
      - name: thioester-His285
        type: angle
        a: {chain: C, resId: 36, atom: O1}
        b: {chain: C, resId: 36, atom: C1}
        c: {chain: A, resId: 285, atom: NE2}

Selections: `chain` plus `resId`/`resIds`, and either `atom`/`atoms` (names), or
`atoms: sidechain | heavy | backbone | all`. Residue numbers follow `chainMap` when one is given,
otherwise the numbering of the topology PDB.

Outputs (in --outDir): geometry_series.csv (one row per frame), geometry_summary.csv / .md
(mean, sd, min, max and the fraction of frames below each --cutoff) and geometry_run.json.
"""
########################################################################################################

BACKBONE = {"N", "CA", "C", "O", "OXT"}


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


def residue_labels(topology, chainMap: Optional[str]) -> Dict[Tuple[str, int], int]:
    """(chainId, residueNumber) -> residue index in the topology."""
    mapping = parse_chain_map(chainMap)
    labels: Dict[Tuple[str, int], int] = {}
    for residue in topology.residues:
        label = None
        for chainId, first, last, offset in mapping:
            if first <= residue.index + 1 <= last:
                label = (chainId, residue.index + 1 - first + 1 + offset)
                break
        if label is None:
            chain = residue.chain
            chainId = getattr(chain, "chain_id", None) or chr(ord("A") + chain.index % 26)
            label = (chainId, residue.resSeq)
        if label in labels:
            raise ValueError(f"residue {label[0]}{label[1]} is not unique; give a --chainMap that separates the chains")
        labels[label] = residue.index
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
        residue = list(topology.residues)[labels[key]]
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


def measure(traj, topology, labels, measurement: dict, massWeighted: bool = True) -> np.ndarray:
    """Return the measurement for every frame, in Angstrom (distances) or degrees (angles)."""
    kind = measurement.get("type", "distance")
    groups = []
    for key in ("a", "b", "c"):
        if key in measurement:
            indices = select_atoms(topology, labels, measurement[key])
            masses = np.array([topology.atom(i).element.mass for i in indices]) if massWeighted else None
            groups.append(centroid(traj.xyz[:, indices, :] * 10.0, masses))
    if kind in ("distance", "comDistance"):
        if len(groups) != 2:
            raise ValueError(f"{measurement.get('name')}: a distance needs selections 'a' and 'b'")
        if kind == "distance":
            for key in ("a", "b"):
                if len(select_atoms(topology, labels, measurement[key])) != 1:
                    raise ValueError(f"{measurement.get('name')}: 'distance' needs one atom per selection; "
                                     f"use type 'comDistance' for a group")
        return np.linalg.norm(groups[0] - groups[1], axis=-1)
    if kind == "angle":
        if len(groups) != 3:
            raise ValueError(f"{measurement.get('name')}: an angle needs selections 'a', 'b' and 'c' (b is the vertex)")
        v1, v2 = groups[0] - groups[1], groups[2] - groups[1]
        cosine = (v1 * v2).sum(-1) / (np.linalg.norm(v1, axis=-1) * np.linalg.norm(v2, axis=-1))
        return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    raise ValueError(f"unknown measurement type '{kind}'")


def summarise(series: pd.DataFrame, cutoffs: Sequence[float]) -> pd.DataFrame:
    rows = []
    for name in series.columns:
        if name in ("frame", "timeNs"):
            continue
        values = series[name].to_numpy()
        row = {"name": name, "mean": values.mean(), "sd": values.std(), "min": values.min(), "max": values.max(),
               "p5": np.percentile(values, 5), "p95": np.percentile(values, 95)}
        for cutoff in cutoffs:
            row[f"frac_lt_{cutoff:g}"] = float((values < cutoff).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def summary_markdown(summary: pd.DataFrame, nFrames: int, cutoffs: Sequence[float]) -> str:
    cutoffColumns = [f"frac_lt_{c:g}" for c in cutoffs]
    header = "| measurement | mean | sd | min | max | " + " | ".join(cutoffColumns) + " |\n"
    header += "|---|---|---|---|---|" + "---|" * len(cutoffColumns) + "\n"
    lines = [f"Geometry over {nFrames} frames (distances in Å, angles in degrees).\n\n", header]
    for _, r in summary.iterrows():
        cells = " | ".join(f"{r[c]:.2f}" for c in cutoffColumns)
        lines.append(f"| {r['name']} | {r['mean']:.2f} | {r['sd']:.2f} | {r['min']:.2f} | {r['max']:.2f} | {cells} |\n")
    return "".join(lines)


########################################################################################################
def main() -> None:
    parser = argparse.ArgumentParser(description="Measure distances, COM distances and angles over a drMD trajectory.")
    parser.add_argument("--pdb", required=True, help="topology PDB")
    parser.add_argument("--trajectory", default=None, help="DCD trajectory; omit to measure the PDB alone")
    parser.add_argument("--measurements", required=True, help="YAML/JSON file describing the measurements")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--frameTimeNs", type=float, default=None, help="time between the used frames (ns)")
    parser.add_argument("--chainMap", default=None, help="overrides the chainMap of the measurement file")
    parser.add_argument("--cutoffs", nargs="*", type=float, default=[3.5, 4.0],
                        help="report the fraction of frames below each of these values")
    parser.add_argument("--geometricCentre", action="store_true", help="use geometric instead of mass-weighted centroids")
    parser.add_argument("--outDir", default="geometry")
    args = parser.parse_args()

    import mdtraj as md

    with open(args.measurements) as fh:
        spec = yaml.safe_load(fh)
    measurements = spec["measurements"] if isinstance(spec, dict) else spec
    chainMap = args.chainMap or (spec.get("chainMap") if isinstance(spec, dict) else None)

    topology = md.load(args.pdb).topology
    traj = md.load(args.trajectory, top=topology, stride=args.stride) if args.trajectory else md.load(args.pdb)
    labels = residue_labels(topology, chainMap)

    dt = args.frameTimeNs
    if dt is None and traj.n_frames > 1:
        dt = float(traj.time[1] - traj.time[0]) / 1000.0
    series = pd.DataFrame({"frame": np.arange(traj.n_frames),
                           "timeNs": np.arange(traj.n_frames) * (dt or 0.0)})
    for measurement in measurements:
        name = measurement["name"]
        series[name] = measure(traj, topology, labels, measurement, massWeighted=not args.geometricCentre)
        print(f"  {name}: mean {series[name].mean():.2f}, range [{series[name].min():.2f}, {series[name].max():.2f}]", flush=True)

    os.makedirs(args.outDir, exist_ok=True)
    series.to_csv(p.join(args.outDir, "geometry_series.csv"), index=False)
    summary = summarise(series, args.cutoffs)
    summary.to_csv(p.join(args.outDir, "geometry_summary.csv"), index=False)
    with open(p.join(args.outDir, "geometry_summary.md"), "w") as fh:
        fh.write(summary_markdown(summary, traj.n_frames, args.cutoffs))
    with open(p.join(args.outDir, "geometry_run.json"), "w") as fh:
        json.dump({"pdb": args.pdb, "trajectory": args.trajectory, "stride": args.stride,
                   "nFrames": int(traj.n_frames), "frameTimeNs": dt, "chainMap": chainMap,
                   "cutoffs": list(args.cutoffs), "massWeighted": not args.geometricCentre,
                   "measurements": measurements}, fh, indent=2)
    print(f"done: {len(measurements)} measurements over {traj.n_frames} frames -> {args.outDir}")


if __name__ == "__main__":
    main()
