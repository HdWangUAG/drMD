"""
Checks for src/ExaminationRoom/drPLIP.py: frame export, ligand atom grouping and summarising.
Run with:  python tests/test_plip.py   (or pytest tests/)
"""
import sys
import tempfile
from os import path as p

import numpy as np
import pandas as pd

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from ExaminationRoom import drPLIP

PDB = p.join(ROOT, "ExampleInputs", "Worked_Example_5_Metadynamics_of_Alanine_Dipeptide", "alanine_dipeptide.pdb")


def test_ppant_group_splits_acyl_chain():
    ## the last three carbons of the acyl chain are the "tail" used as a metadynamics CV group
    assert drPLIP.ppant_group("C12", "S12") == "tail"
    assert drPLIP.ppant_group("C10", "S12") == "tail"
    assert drPLIP.ppant_group("C9", "S12") == "acyl"
    assert drPLIP.ppant_group("C12", "S14") == "tail"
    assert drPLIP.ppant_group("C11", "S14") == "acyl"
    ## named groups win over the carbon-number rule
    assert drPLIP.ppant_group("C1", "S12") == "thioester"
    assert drPLIP.ppant_group("O1", "S12") == "thioester"
    assert drPLIP.ppant_group("P24", "S12") == "phosphate"
    assert drPLIP.ppant_group("C32", "S12") == "pant"
    assert drPLIP.ppant_group("OG", "S12") == "ser"


def test_graph_groups_classifies_a_thioester_without_atom_names():
    ## S - C1(=O1) - C2 - ... - C6 : a six-carbon acyl chain, arbitrary atom names
    coords = {1: ("X1", "LIG", "E", 1, "S", (0.0, 0.0, 0.0)),
              2: ("X2", "LIG", "E", 1, "C", (1.8, 0.0, 0.0)),
              3: ("X3", "LIG", "E", 1, "O", (2.4, 1.1, 0.0)),
              4: ("X4", "LIG", "E", 1, "C", (3.3, 0.0, 0.0)),
              5: ("X5", "LIG", "E", 1, "C", (4.8, 0.0, 0.0)),
              6: ("X6", "LIG", "E", 1, "C", (6.3, 0.0, 0.0)),
              7: ("X7", "LIG", "E", 1, "C", (7.8, 0.0, 0.0)),
              8: ("X8", "LIG", "E", 1, "C", (9.3, 0.0, 0.0)),
              9: ("X9", "LIG", "E", 1, "P", (-3.0, 0.0, 0.0)),
              10: ("XA", "LIG", "E", 1, "O", (-4.5, 0.0, 0.0))}
    groups = drPLIP.graph_groups(coords, list(coords))
    assert groups[1] == groups[2] == groups[3] == "thioester"
    assert groups[4] == groups[5] == "acyl"      # C2, C3 of the six-carbon chain
    assert groups[6] == groups[7] == groups[8] == "tail"   # last three carbons C4, C5, C6
    assert groups[9] == groups[10] == "phosphate"


def test_parse_chain_map():
    entries = drPLIP.parse_chain_map("A:1-299:83,C:599-675:0")
    assert entries == [("A", 1, 299, 83), ("C", 599, 675, 0)]
    assert drPLIP.parse_chain_map(None) == []


def test_write_frames_strips_hydrogens_and_renumbers():
    with tempfile.TemporaryDirectory() as tmp:
        frames = drPLIP.write_frames(PDB, None, tmp, ligandResidues=["NME"], chainMap="A:1-3:100")
        assert len(frames) == 1
        atoms = drPLIP.read_pdb_atoms(frames[0][0])
        assert atoms, "no atoms written"
        assert all(a[4] != "H" for a in atoms.values()), "hydrogens were not stripped"
        ## residue 1 of the topology becomes chain A, number 101 (1 + offset 100)
        assert {a[2] for a in atoms.values()} == {"A"}
        assert min(a[3] for a in atoms.values()) == 101
        ## the residue named as a ligand is written as HETATM
        text = open(frames[0][0]).read()
        assert any(l.startswith("HETATM") and l[17:20] == "NME" for l in text.splitlines())
        assert any(l.startswith("ATOM") and l[17:20] == "ALA" for l in text.splitlines())


def test_nearest_atom_finds_the_reported_coordinate():
    xyz = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    serials = [11, 22, 33]
    assert drPLIP.nearest_atom(xyz, serials, np.array([4.9, 0.1, 0.0])) == 22
    assert drPLIP.nearest_atom(xyz, serials, None) == -1


def test_to_float_handles_numpy_reprs():
    assert drPLIP._to_float("np.float64(43.385)") == 43.385
    assert drPLIP._to_float("3.71") == 3.71
    coords = drPLIP._parse_coords("(np.float64(1.0), np.float64(-2.5), np.float64(3.0))")
    assert list(coords) == [1.0, -2.5, 3.0]


def test_summarise_counts_frames_not_rows():
    ## the same contact seen through two atoms in one frame must not count twice
    rows = []
    for frame in range(4):
        n = 2 if frame < 2 else 1
        for _ in range(n):
            rows.append(dict(view="ligand", site="S12:C:36", itype="hydrophobic", partnerChain="A",
                             partnerResNum=189, partnerResName="ILE", ligandChain="C", ligandResNum=36,
                             ligandResName="S12", ligandGroup="acyl", frame=frame, dist=3.6))
    rows.append(dict(view="ligand", site="S12:C:36", itype="hbond", partnerChain="A", partnerResNum=285,
                     partnerResName="HIS", ligandChain="C", ligandResNum=36, ligandResName="S12",
                     ligandGroup="thioester", frame=0, dist=2.9))
    summary = drPLIP.summarise(pd.DataFrame(rows), nFrames=4)
    hydrophobic = summary[summary.itype == "hydrophobic"].iloc[0]
    assert hydrophobic.nFrames == 4 and hydrophobic.fraction == 1.0
    hbond = summary[summary.itype == "hbond"].iloc[0]
    assert hbond.fraction == 0.25
    ## minFraction filters the rare one out
    assert drPLIP.summarise(pd.DataFrame(rows), nFrames=4, minFraction=0.5).itype.tolist() == ["hydrophobic"]
    assert "HIS 285:A" in drPLIP.summary_markdown(summary, 4, 0.0)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all drPLIP tests passed")
