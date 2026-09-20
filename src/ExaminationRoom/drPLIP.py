## BASIC PYTHON LIBRARIES
import os
from os import path as p
import argparse
import csv
import json
import re

## NUMERICAL LIBRARIES
import numpy as np
import pandas as pd

## CLEAN CODE
from typing import Dict, List, Optional, Sequence, Tuple

########################################################################################################
"""
Interaction profiling of drMD trajectories with PLIP (Protein-Ligand Interaction Profiler).

Two views are produced for every frame:

  ligand view   : each residue named in --ligandResidues (e.g. an acyl-Ppant-Ser ncAA, or a
                  HETATM ligand) is treated as a ligand; PLIP reports its non-covalent contacts
                  with every other residue (KEEPMOD mode, so modified residues are kept as ligands).
  peptide view  : each chain named in --peptideChains is treated as a peptide ligand; PLIP reports
                  the residue-residue contacts across the protein-protein interface.

PLIP re-protonates the structure with Open Babel, so hydrogens are stripped from the frames before
profiling: hydrogen-bond assignments therefore reflect Open Babel's protonation, NOT the tautomers
of the MD topology. Water is not written, so PLIP water bridges are not available here.

Outputs (in --outDir):
  plip_interactions.csv   one row per interaction per frame (long format)
  plip_summary.csv        per (partner residue, ligand residue/group, interaction type): fraction of
                          frames in which the interaction is present, mean distance
  plip_summary.md         the same, as tables for a report
  frames/                 the PDB files handed to PLIP (kept for inspection / PyMOL)

PLIP: Adasme et al., Nucleic Acids Res. 2021; Schake et al., Nucleic Acids Res. 2025 (PLIP 2025).
"""
########################################################################################################

PLIP_TABLES: List[str] = ["hydrophobic", "hbond", "waterbridge", "saltbridge",
                          "pistacking", "pication", "halogen", "metal"]

## acyl-Ppant-Ser atom groups (S12/S14 templates of the MD_HFA campaign); other ligands fall back to "ligand"
PPANT_GROUPS: Dict[str, set] = {
    "ser": {"N", "CA", "C", "O", "CB", "OG"},
    "phosphate": {"P24", "O23", "O26", "O27"},
    "pant": {"C28", "C29", "C30", "C31", "C32", "O33", "C34", "O35", "N36",
             "C37", "C38", "C39", "O40", "N41", "C42", "C43"},
    "thioester": {"S1", "C1", "O1"},
}


## a number that is not part of an identifier: "np.float64(43.385)" must give 43.385, not 64
NUMBER = re.compile(r"(?<![A-Za-z0-9_.])-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def _to_float(value) -> float:
    """PLIP report cells may be plain numbers or repr'd numpy scalars."""
    m = NUMBER.search(str(value))
    return float(m.group()) if m else float("nan")


def _parse_coords(value) -> Optional[np.ndarray]:
    nums = NUMBER.findall(str(value))
    if len(nums) < 3:
        return None
    return np.array([float(x) for x in nums[:3]])


def ppant_group(atomName: str, resName: str) -> str:
    """Classify an atom of an S<n> acyl-Ppant-Ser residue: ser / phosphate / pant / thioester / acyl / tail."""
    for group, names in PPANT_GROUPS.items():
        if atomName in names:
            return group
    m = re.fullmatch(r"S(\d+)", resName)
    if m and re.fullmatch(r"C\d+", atomName):
        nCarbons = int(m.group(1))
        k = int(atomName[1:])
        return "tail" if k > nCarbons - 3 else "acyl"
    return "ligand"


def graph_groups(atoms: Dict[int, tuple], serialsOfLigand: List[int]) -> Dict[int, str]:
    """Classify a thioester-acyl-Ppant ligand without an atom-name convention (e.g. Boltz output):
    phosphate = P and its O; thioester = S, the carbonyl C bonded to it and its O;
    acyl/tail = the carbon chain walked from the carbonyl carbon away from S; everything else = pant."""
    xyz = np.array([atoms[s][5] for s in serialsOfLigand])
    dist = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    neighbours = {s: [serialsOfLigand[j] for j in range(len(serialsOfLigand)) if i != j and dist[i, j] < 1.9]
                  for i, s in enumerate(serialsOfLigand)}
    element = lambda s: atoms[s][4]
    groups: Dict[int, str] = {}
    for s in serialsOfLigand:
        if element(s) == "P":
            groups[s] = "phosphate"
            for o in neighbours[s]:
                groups[o] = "phosphate"
    for s in serialsOfLigand:
        if element(s) != "S":
            continue
        groups[s] = "thioester"
        carbonyls = [c for c in neighbours[s] if element(c) == "C" and any(element(o) == "O" for o in neighbours[c])]
        if not carbonyls:
            continue
        c1 = carbonyls[0]
        groups[c1] = "thioester"
        for o in neighbours[c1]:
            if element(o) == "O":
                groups[o] = "thioester"
        chain, prev, cur = [c1], s, c1
        while True:
            nxt = [c for c in neighbours[cur] if element(c) == "C" and c != prev and c not in chain]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            chain.append(cur)
        for k, c in enumerate(chain[1:], start=2):
            groups[c] = "tail" if k > len(chain) - 3 else "acyl"
    for s in serialsOfLigand:
        groups.setdefault(s, "pant")
    return groups


########################################################################################################
def parse_chain_map(chainMap: Optional[str]) -> List[Tuple[str, int, int, int]]:
    """'A:1-299:83,B:300-598:83,C:599-675:0' -> [(chainId, firstResidueIndex1, lastResidueIndex1, numberingOffset)]
    Residue indices are 1-based positions in the topology; the offset is added to the within-chain
    position (1-based) to give the residue number written to the PDB."""
    if not chainMap:
        return []
    entries = []
    for item in chainMap.split(","):
        chainId, span, offset = item.strip().split(":")
        first, last = span.split("-")
        entries.append((chainId, int(first), int(last), int(offset)))
    return entries


def write_frames(pdbFile: str,
                 trajectory: Optional[str],
                 outDir: str,
                 ligandResidues: Sequence[str],
                 chainMap: Optional[str] = None,
                 stride: int = 1,
                 frameTimeNs: Optional[float] = None,
                 keepSolvent: bool = False) -> List[Tuple[str, float]]:
    """Write hydrogen-free PDB frames for PLIP. Returns [(pdbPath, timeNs)]."""
    import mdtraj as md

    top = md.load(pdbFile).topology
    if trajectory:
        traj = md.load(trajectory, top=top, stride=stride)
    else:
        traj = md.load(pdbFile)
    if not keepSolvent:
        keep = [a.index for a in top.atoms if not a.residue.is_water and a.residue.name not in
                ("NA", "CL", "K", "Na+", "Cl-", "K+", "MG", "ZN", "CA")]
        traj = traj.atom_slice(keep)
        top = traj.topology

    mapping = parse_chain_map(chainMap)
    chainOf: Dict[int, Tuple[str, int]] = {}
    for residue in top.residues:
        assigned = False
        for chainId, first, last, offset in mapping:
            if first <= residue.index + 1 <= last:
                chainOf[residue.index] = (chainId, residue.index + 1 - first + 1 + offset)
                assigned = True
                break
        if not assigned:
            chainId = residue.chain.chain_id if getattr(residue.chain, "chain_id", None) else chr(ord("A") + residue.chain.index % 26)
            chainOf[residue.index] = (chainId, residue.resSeq)

    os.makedirs(outDir, exist_ok=True)
    fmt = ("{rec:<6}{serial:>5} {name:<4}{alt:1}{res:>3} {ch:1}{seq:>4}{ic:1}   "
           "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {el:>2}\n")
    stem = p.splitext(p.basename(trajectory or pdbFile))[0]
    written = []
    dt = frameTimeNs if frameTimeNs is not None else (traj.time[1] - traj.time[0]) / 1000.0 if traj.n_frames > 1 else 0.0
    for f in range(traj.n_frames):
        xyz = traj.xyz[f] * 10.0
        lines, serial, prevChain = [], 0, None
        for atom in top.atoms:
            if atom.element is not None and atom.element.symbol == "H":
                continue
            residue = atom.residue
            chainId, resSeq = chainOf[residue.index]
            if prevChain is not None and chainId != prevChain:
                lines.append("TER\n")
            prevChain = chainId
            serial += 1
            name = atom.name if len(atom.name) == 4 else " " + atom.name
            lines.append(fmt.format(rec="HETATM" if residue.name in ligandResidues else "ATOM",
                                    serial=serial % 100000, name=name, alt="", res=residue.name[:3],
                                    ch=chainId, seq=resSeq, ic="",
                                    x=xyz[atom.index, 0], y=xyz[atom.index, 1], z=xyz[atom.index, 2],
                                    occ=1.0, b=0.0, el=atom.element.symbol if atom.element else name.strip()[0]))
        lines.append("TER\nEND\n")
        timeNs = f * dt * (1 if trajectory else 0)
        outPdb = p.join(outDir, f"{stem}_f{f:04d}.pdb")
        with open(outPdb, "w") as fh:
            fh.writelines(lines)
        written.append((outPdb, timeNs))
    return written


########################################################################################################
def read_pdb_atoms(pdbFile: str) -> Dict[int, tuple]:
    """serial -> (atomName, resName, chainId, resNum, element, xyz)"""
    atoms = {}
    for line in open(pdbFile):
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            element = line[76:78].strip() or re.sub(r"[^A-Za-z]", "", name)[:1]
            atoms[int(line[6:11])] = (name, line[17:20].strip(), line[21], int(line[22:26]), element.upper(),
                                      (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def run_plip(pdbFile: str, peptideChains: Sequence[str], keepModified: bool) -> List[Tuple[str, str, dict]]:
    """Run PLIP once; returns [(bindingSiteKey, table, rowDict)]."""
    from plip.basic import config
    from plip.structure.preparation import PDBComplex
    from plip.exchange.report import BindingSiteReport

    config.PEPTIDES = list(peptideChains)
    config.KEEPMOD = keepModified
    config.NOHYDRO = False
    config.INTRA = None
    cx = PDBComplex()
    cx.load_pdb(pdbFile)
    for ligand in cx.ligands:
        cx.characterize_complex(ligand)
    rows = []
    for key, site in cx.interaction_sets.items():
        report = BindingSiteReport(site)
        for table in PLIP_TABLES:
            features = getattr(report, f"{table}_features")
            for info in getattr(report, f"{table}_info"):
                rows.append((key, table, dict(zip(features, info))))
    return rows


def nearest_atom(atomsXyz: np.ndarray, serials: List[int], coords: Optional[np.ndarray]) -> int:
    """PLIP's atom indices are renumbered internally (TER records count), so ligand atoms are
    identified from the reported coordinates instead. Charge-group / ring centroids map to the
    nearest member atom."""
    if coords is None:
        return -1
    return serials[int(np.argmin(((atomsXyz - coords) ** 2).sum(1)))]


INTERACTION_FIELDS = ["label", "frame", "timeNs", "view", "site", "itype",
                      "partnerChain", "partnerResNum", "partnerResName",
                      "ligandChain", "ligandResNum", "ligandResName", "ligandAtom", "ligandGroup",
                      "dist", "distDA", "angle", "partnerSidechain", "partnerIsDonor", "partnerIsPositive"]


def profile_one_frame(job: Tuple[int, str, float, str, Sequence[str], Sequence[str]]) -> List[dict]:
    """One frame, both views. A free function so that it can run in a worker process."""
    frameIndex, pdbFile, timeNs, label, ligandResidues, peptideChains = job
    atoms = read_pdb_atoms(pdbFile)
    serials = list(atoms)
    atomsXyz = np.array([atoms[s][5] for s in serials])
    graphGroups: Dict[str, Dict[int, str]] = {}
    views = [("ligand", [], True)]
    if peptideChains:
        views.append(("peptide", list(peptideChains), False))
    rows = []
    for view, chains, keepMod in views:
        for site, table, d in run_plip(pdbFile, chains, keepMod):
            serial = nearest_atom(atomsXyz, serials, _parse_coords(d.get("LIGCOO")))
            atomName, group = "", ""
            if serial in atoms:
                atomName, resName, ligChain = atoms[serial][0], atoms[serial][1], atoms[serial][2]
                if resName in ligandResidues and re.fullmatch(r"S\d+", resName):
                    group = ppant_group(atomName, resName)
                elif resName in ligandResidues or d["RESTYPE_LIG"] in ligandResidues:
                    key = f"{ligChain}:{atoms[serial][3]}"
                    if key not in graphGroups:
                        members = [s for s in serials if atoms[s][2] == ligChain and atoms[s][3] == atoms[serial][3]]
                        graphGroups[key] = graph_groups(atoms, members)
                    group = graphGroups[key].get(serial, "ligand")
                else:
                    group = "protein"
            dist = d.get("DIST", d.get("DIST_H-A", d.get("CENTDIST", "")))
            rows.append(dict(zip(INTERACTION_FIELDS, [
                label, frameIndex, timeNs, view, site, table,
                d["RESCHAIN"], int(d["RESNR"]), d["RESTYPE"],
                d["RESCHAIN_LIG"], int(d["RESNR_LIG"]), d["RESTYPE_LIG"], atomName, group,
                _to_float(dist), _to_float(d.get("DIST_D-A", "nan")),
                _to_float(d.get("DON_ANGLE", d.get("ANGLE", "nan"))),
                d.get("SIDECHAIN", ""), d.get("PROTISDON", ""), d.get("PROTISPOS", "")])))
    return rows


def profile_frames(frames: Sequence[Tuple[str, float]],
                   label: str,
                   ligandResidues: Sequence[str],
                   peptideChains: Sequence[str],
                   outCsv: str,
                   nProc: int = 1) -> pd.DataFrame:
    """PLIP is single-threaded, so frames are profiled in parallel processes."""
    jobs = [(i, pdbFile, timeNs, label, list(ligandResidues), list(peptideChains))
            for i, (pdbFile, timeNs) in enumerate(frames)]
    rows: List[dict] = []
    if nProc > 1 and len(jobs) > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(min(nProc, len(jobs))) as pool:
            for done, frameRows in enumerate(pool.imap_unordered(profile_one_frame, jobs, chunksize=1), start=1):
                rows.extend(frameRows)
                print(f"  PLIP {done}/{len(jobs)}", flush=True)
    else:
        for job in jobs:
            rows.extend(profile_one_frame(job))
            print(f"  PLIP {job[0] + 1}/{len(jobs)} {p.basename(job[1])}", flush=True)
    rows.sort(key=lambda r: (r["frame"], r["view"], r["site"], r["itype"]))
    df = pd.DataFrame(rows, columns=INTERACTION_FIELDS)
    df.to_csv(outCsv, index=False)
    return df


########################################################################################################
def summarise(df: pd.DataFrame, nFrames: int, minFraction: float = 0.0) -> pd.DataFrame:
    """Fraction of frames in which each (partner residue, ligand residue/group, type) occurs."""
    if df.empty:
        return pd.DataFrame()
    key = ["view", "site", "itype", "partnerChain", "partnerResNum", "partnerResName",
           "ligandChain", "ligandResNum", "ligandResName", "ligandGroup"]
    perFrame = df.drop_duplicates(key + ["frame"])
    counts = perFrame.groupby(key).agg(nFrames=("frame", "nunique"),
                                       meanDist=("dist", "mean"),
                                       minDist=("dist", "min")).reset_index()
    counts["fraction"] = counts["nFrames"] / nFrames
    counts = counts[counts["fraction"] >= minFraction]
    return counts.sort_values(["view", "site", "fraction"], ascending=[True, True, False])


def summary_markdown(summary: pd.DataFrame, nFrames: int, minFraction: float) -> str:
    out = [f"PLIP interaction frequencies over {nFrames} frames (interactions present in >= {minFraction:.0%} of frames).\n"]
    if summary.empty:
        return "".join(out) + "\nno interactions found\n"
    for (view, site), block in summary.groupby(["view", "site"], sort=False):
        out.append(f"\n### {view} view — {site}\n\n")
        out.append("| type | partner | ligand residue / group | fraction | mean dist (Å) |\n|---|---|---|---|---|\n")
        for _, r in block.iterrows():
            partner = f"{r.partnerResName} {r.partnerResNum}:{r.partnerChain}"
            lig = f"{r.ligandResName} {r.ligandResNum}:{r.ligandChain}"
            if r.ligandGroup and r.ligandGroup != "protein":
                lig += f" [{r.ligandGroup}]"
            out.append(f"| {r.itype} | {partner} | {lig} | {r.fraction:.2f} | {r.meanDist:.2f} |\n")
    return "".join(out)


########################################################################################################
def main() -> None:
    parser = argparse.ArgumentParser(description="Profile protein-ligand and protein-protein interactions of a drMD trajectory with PLIP.")
    parser.add_argument("--pdb", required=True, help="topology PDB (e.g. <step>/trajectory.pdb or a solute-only PDB)")
    parser.add_argument("--trajectory", default=None, help="DCD trajectory; omit to profile the PDB alone")
    parser.add_argument("--stride", type=int, default=1, help="use every n-th frame")
    parser.add_argument("--frameTimeNs", type=float, default=None, help="time between written frames (ns), for the output tables")
    parser.add_argument("--ligandResidues", nargs="*", default=[], help="residue names to treat as ligands (written as HETATM)")
    parser.add_argument("--peptideChains", nargs="*", default=[], help="chain IDs to treat as peptide ligands (protein-protein view)")
    parser.add_argument("--chainMap", default=None,
                        help="restore chain IDs / numbering from residue positions, e.g. 'A:1-299:83,B:300-598:83,C:599-675:0'")
    parser.add_argument("--label", default="run", help="label column in the output tables")
    parser.add_argument("--minFraction", type=float, default=0.1, help="report interactions present in at least this fraction of frames")
    parser.add_argument("--outDir", default="plip")
    parser.add_argument("--nProc", type=int, default=1, help="profile this many frames in parallel (PLIP itself is single-threaded)")
    parser.add_argument("--keepSolvent", action="store_true", help="keep water and ions in the frames")
    args = parser.parse_args()

    os.makedirs(args.outDir, exist_ok=True)
    frames = write_frames(args.pdb, args.trajectory, p.join(args.outDir, "frames"), args.ligandResidues,
                          chainMap=args.chainMap, stride=args.stride, frameTimeNs=args.frameTimeNs,
                          keepSolvent=args.keepSolvent)
    print(f"wrote {len(frames)} frames")
    df = profile_frames(frames, args.label, args.ligandResidues, args.peptideChains,
                        p.join(args.outDir, "plip_interactions.csv"), nProc=args.nProc)
    summary = summarise(df, len(frames), args.minFraction)
    summary.to_csv(p.join(args.outDir, "plip_summary.csv"), index=False)
    with open(p.join(args.outDir, "plip_summary.md"), "w") as fh:
        fh.write(summary_markdown(summary, len(frames), args.minFraction))
    with open(p.join(args.outDir, "plip_run.json"), "w") as fh:
        json.dump({"pdb": args.pdb, "trajectory": args.trajectory, "stride": args.stride, "nFrames": len(frames),
                   "ligandResidues": list(args.ligandResidues), "peptideChains": list(args.peptideChains),
                   "chainMap": args.chainMap, "nProc": args.nProc, "note": "hydrogens stripped; PLIP/Open Babel protonation; no water"},
                  fh, indent=2)
    print(f"done: {len(df)} interactions, summary in {args.outDir}/plip_summary.md")


if __name__ == "__main__":
    main()
