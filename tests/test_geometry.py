"""
Checks for src/ExaminationRoom/drGeometry.py.
Run with:  python tests/test_geometry.py   (or pytest tests/)
"""
import sys
from os import path as p

import numpy as np

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from ExaminationRoom import drGeometry

import mdtraj as md

PDB = p.join(ROOT, "ExampleInputs", "Worked_Example_5_Metadynamics_of_Alanine_Dipeptide", "alanine_dipeptide.pdb")

TRAJ = md.load(PDB)
TOP = TRAJ.topology
LABELS = drGeometry.residue_labels(TOP, None)
RESIDUES = [(chainId, number) for chainId, number in LABELS]


def test_residue_labels_are_unique_and_cover_the_topology():
    assert len(LABELS) == TOP.n_residues
    ## a chainMap renumbers: residue 1 of the topology becomes 101
    labels = drGeometry.residue_labels(TOP, f"A:1-{TOP.n_residues}:100")
    assert ("A", 101) in labels and labels[("A", 101)] == 0


def test_select_atoms_by_name_and_by_class():
    chainId, number = RESIDUES[1]
    one = drGeometry.select_atoms(TOP, LABELS, {"chain": chainId, "resId": number, "atom": "CA"})
    assert len(one) == 1 and TOP.atom(one[0]).name == "CA"
    heavy = drGeometry.select_atoms(TOP, LABELS, {"chain": chainId, "resId": number, "atoms": "heavy"})
    assert all(TOP.atom(i).element.symbol != "H" for i in heavy)
    backbone = drGeometry.select_atoms(TOP, LABELS, {"chain": chainId, "resId": number, "atoms": "backbone"})
    assert {TOP.atom(i).name for i in backbone} <= drGeometry.BACKBONE
    sidechain = drGeometry.select_atoms(TOP, LABELS, {"chain": chainId, "resId": number, "atoms": "sidechain"})
    assert not ({TOP.atom(i).name for i in sidechain} & drGeometry.BACKBONE)
    assert set(backbone) | set(sidechain) == set(heavy)
    named = drGeometry.select_atoms(TOP, LABELS, {"chain": chainId, "resId": number, "atoms": ["N", "CA"]})
    assert len(named) == 2


def test_select_atoms_reports_missing_atoms_and_residues():
    chainId, number = RESIDUES[0]
    for selection, error in [({"chain": chainId, "resId": number, "atom": "ZZ"}, ValueError),
                             ({"chain": chainId, "resId": 9999, "atom": "CA"}, KeyError),
                             ({"resId": number, "atom": "CA"}, ValueError)]:
        try:
            drGeometry.select_atoms(TOP, LABELS, selection)
        except error:
            continue
        raise AssertionError(f"expected {error.__name__} for {selection}")


def test_distance_matches_mdtraj():
    (chain1, res1), (chain2, res2) = RESIDUES[0], RESIDUES[2]
    a = drGeometry.select_atoms(TOP, LABELS, {"chain": chain1, "resId": res1, "atom": "C"})[0]
    b = drGeometry.select_atoms(TOP, LABELS, {"chain": chain2, "resId": res2, "atom": "N"})[0]
    reference = md.compute_distances(TRAJ, [[a, b]], periodic=False)[0, 0] * 10.0
    value = drGeometry.measure(TRAJ, TOP, LABELS, {"name": "d", "type": "distance",
                                                   "a": {"chain": chain1, "resId": res1, "atom": "C"},
                                                   "b": {"chain": chain2, "resId": res2, "atom": "N"}})
    assert np.isclose(value[0], reference, atol=1e-4)


def test_distance_rejects_multi_atom_selections():
    chainId, number = RESIDUES[1]
    spec = {"name": "d", "type": "distance",
            "a": {"chain": chainId, "resId": number, "atoms": "heavy"},
            "b": {"chain": chainId, "resId": number, "atom": "CA"}}
    try:
        drGeometry.measure(TRAJ, TOP, LABELS, spec)
    except ValueError as error:
        assert "comDistance" in str(error)
        return
    raise AssertionError("a multi-atom 'distance' should be refused")


def test_com_distance_is_mass_weighted_by_default():
    chainId, number = RESIDUES[1]
    other, otherNumber = RESIDUES[2]
    spec = {"name": "com", "type": "comDistance",
            "a": {"chain": chainId, "resId": number, "atoms": "sidechain"},
            "b": {"chain": other, "resId": otherNumber, "atoms": "heavy"}}
    massWeighted = drGeometry.measure(TRAJ, TOP, LABELS, spec)[0]
    geometric = drGeometry.measure(TRAJ, TOP, LABELS, spec, massWeighted=False)[0]
    assert massWeighted > 0 and geometric > 0
    ## the two centroids differ, so the distances must not be identical
    assert not np.isclose(massWeighted, geometric, atol=1e-6)


def test_angle_of_three_atoms():
    chainId, number = RESIDUES[1]
    spec = {"name": "angle", "type": "angle",
            "a": {"chain": chainId, "resId": number, "atom": "N"},
            "b": {"chain": chainId, "resId": number, "atom": "CA"},
            "c": {"chain": chainId, "resId": number, "atom": "C"}}
    value = drGeometry.measure(TRAJ, TOP, LABELS, spec)[0]
    indices = [drGeometry.select_atoms(TOP, LABELS, spec[k])[0] for k in ("a", "b", "c")]
    reference = np.degrees(md.compute_angles(TRAJ, [indices], periodic=False)[0, 0])
    assert np.isclose(value, reference, atol=1e-3)
    assert 80.0 < value < 140.0


def test_summarise_reports_cutoff_fractions():
    import pandas as pd
    series = pd.DataFrame({"frame": range(4), "timeNs": [0.0, 0.2, 0.4, 0.6], "d": [2.0, 3.0, 4.0, 5.0]})
    summary = drGeometry.summarise(series, [3.5, 4.5])
    row = summary.iloc[0]
    assert row["name"] == "d" and row["mean"] == 3.5 and row["min"] == 2.0 and row["max"] == 5.0
    assert row["frac_lt_3.5"] == 0.5 and row["frac_lt_4.5"] == 0.75
    assert "frac_lt_3.5" in drGeometry.summary_markdown(summary, 4, [3.5, 4.5])


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all drGeometry tests passed")
