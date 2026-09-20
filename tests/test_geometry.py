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



########################################################################################################
## a synthetic solvated system, so that every expected answer can be worked out by hand
from mdtraj.core import element

BOX = 30.0  ## Angstrom, cubic: big enough that only the water we place across a face wraps


def build_system(waterOxygens, box=BOX, waterName="HOH"):
    """Two solute residues (S12 C1/O1, HIS NE2) plus waters. waterOxygens: (frames, nWaters, 3) in Angstrom."""
    topology = md.Topology()
    ligChain = topology.add_chain()
    ligand = topology.add_residue("S12", ligChain, resSeq=1)
    topology.add_atom("C1", element.carbon, ligand)
    topology.add_atom("O1", element.oxygen, ligand)
    hisChain = topology.add_chain()
    his = topology.add_residue("HIS", hisChain, resSeq=1)
    topology.add_atom("NE2", element.nitrogen, his)
    waterChain = topology.add_chain()
    waterOxygens = np.asarray(waterOxygens, dtype=float)
    nFrames, nWaters, _ = waterOxygens.shape
    for i in range(nWaters):
        water = topology.add_residue(waterName, waterChain, resSeq=i + 1)
        topology.add_atom("O", element.oxygen, water)
        topology.add_atom("H1", element.hydrogen, water)
        topology.add_atom("H2", element.hydrogen, water)

    solute = np.array([[0.0, 0.0, 0.0], [1.23, 0.0, 0.0], [6.0, 0.0, 0.0]])
    xyz = np.zeros((nFrames, topology.n_atoms, 3))
    xyz[:, 0:3, :] = solute[None, :, :]
    for i in range(nWaters):
        base = 3 + 3 * i
        xyz[:, base, :] = waterOxygens[:, i, :]
        ## hydrogens are never measured, they only have to exist
        xyz[:, base + 1, :] = waterOxygens[:, i, :] + np.array([0.96, 0.0, 0.0])
        xyz[:, base + 2, :] = waterOxygens[:, i, :] + np.array([-0.24, 0.93, 0.0])
    traj = md.Trajectory(xyz=xyz / 10.0, topology=topology,
                         unitcell_lengths=np.tile([box / 10.0] * 3, (nFrames, 1)),
                         unitcell_angles=np.tile([90.0] * 3, (nFrames, 1)))
    return traj


##                       W1 (nearest in frame 0)  W2 (bridges, nearest in frame 1)  W3 (near C1 only)
##                       W4 (bulk)                W5 (bridges in frame 1 only)      W6 (near NE2 only)
WATERS = np.array([[[0.0, 3.0, 0.0], [3.0, 1.0, 0.0], [0.0, 3.9, 0.0],
                    [10.0, 10.0, 10.0], [0.0, 0.0, 12.0], [6.0, 2.0, 0.0]],
                   [[0.0, 12.0, 0.0], [3.0, 1.0, 0.0], [0.0, 3.9, 0.0],
                    [10.0, 10.0, 10.0], [3.2, -1.0, 0.0], [6.0, 2.0, 0.0]]])
SOLVENT_TRAJ = build_system(WATERS)
SOLVENT_TOP = SOLVENT_TRAJ.topology
SOLVENT_LABELS = drGeometry.residue_labels(SOLVENT_TOP, None)
C1 = {"chain": "A", "resId": 1, "atom": "C1"}
O1 = {"chain": "A", "resId": 1, "atom": "O1"}
NE2 = {"chain": "B", "resId": 1, "atom": "NE2"}


def solvent_measure(spec, traj=None):
    return drGeometry.measure(traj if traj is not None else SOLVENT_TRAJ, SOLVENT_TOP, SOLVENT_LABELS, spec)


def test_solvent_oxygens_are_found_by_is_water_and_by_name():
    oxygens = drGeometry.solvent_oxygen_indices(SOLVENT_TOP)
    assert list(oxygens) == [3, 6, 9, 12, 15, 18]
    ## a tleap-style name that mdtraj's is_water does not know
    t3p = build_system(WATERS, waterName="T3P")
    assert not list(t3p.topology.residues)[2].is_water
    assert drGeometry.solvent_oxygen_indices(t3p.topology).size == 6
    ## and the set is configurable, so a name can be excluded too
    assert drGeometry.solvent_oxygen_indices(t3p.topology, ["HOH"]).size == 0


def test_nearest_solvent_picks_the_right_water_in_each_frame():
    values = solvent_measure({"name": "nearest", "type": "nearestSolvent", "a": C1})
    ## frame 0: W1 at 3.0 A; frame 1: W1 has left and W2 at sqrt(10) A is closest
    assert np.allclose(values, [3.0, np.sqrt(10.0)], atol=1e-4)
    ## a multi-atom selection measures from its centroid: C1/O1 midpoint is 0.615 A along x
    both = solvent_measure({"name": "nearestGroup", "type": "nearestSolvent",
                            "a": {"chain": "A", "resId": 1, "atoms": ["C1", "O1"]}})
    assert both[0] > 0 and both[0] < 4.0


def test_solvent_bridge_counts_only_waters_satisfying_both_cutoffs():
    spec = {"name": "bridge", "type": "solventBridge", "a": C1, "b": NE2, "cutoffA": 4.0, "cutoffB": 3.5}
    counts = solvent_measure(spec)
    ## frame 0: only W2 (3.16 A of C1, 3.16 A of NE2). W3 is within 4 A of C1 but 7.2 A of NE2,
    ## W6 is within 3.5 A of NE2 but 6.3 A of C1; frame 1 adds W5.
    assert list(counts) == [1.0, 2.0]
    ## tightening the first cutoff below sqrt(10) must drop W2 and leave nothing in frame 0
    tight = solvent_measure({**spec, "cutoffA": 3.0})
    assert list(tight) == [0.0, 0.0]
    ## the defaults are 4.0 / 3.5
    assert list(solvent_measure({"name": "bridge", "type": "solventBridge", "a": C1, "b": NE2})) == [1.0, 2.0]


def test_bridge_angle_matches_a_hand_computed_angle():
    ## Burgi-Dunitz style: Ow ... C1=O1, with the carbonyl carbon as the vertex
    values = solvent_measure({"name": "bd", "type": "bridgeAngle", "a": C1, "b": C1, "c": O1})
    ## frame 0 the nearest water sits on +y, exactly perpendicular to C1->O1 (+x); frame 1 it is at
    ## (3, 1, 0), so the angle is arccos(3 / sqrt(10))
    expected = [90.0, np.degrees(np.arccos(3.0 / np.sqrt(10.0)))]
    assert np.allclose(values, expected, atol=1e-3)


def test_per_frame_selection_follows_the_frame_not_the_first_one():
    """The diagonal test: a different water is nearest in each frame, so a [:, 0] bug would show up."""
    values = solvent_measure({"name": "reach", "type": "nearestSolventTo", "a": C1, "b": NE2})
    ## frame 0 follows W1 at (0, 3, 0) -> sqrt(45) from NE2; frame 1 follows W2 at (3, 1, 0) -> sqrt(10)
    expected = [np.sqrt(45.0), np.sqrt(10.0)]
    assert np.allclose(values, expected, atol=1e-4)
    ## the two frames must genuinely differ, otherwise this test would pass with the bug in place
    assert abs(expected[0] - expected[1]) > 3.0
    ## and each frame must equal the same measurement made on that frame alone
    for frame in (0, 1):
        single = solvent_measure({"name": "reach", "type": "nearestSolventTo", "a": C1, "b": NE2},
                                 traj=SOLVENT_TRAJ[frame])
        assert np.isclose(single[0], values[frame], atol=1e-6)


def test_per_frame_distances_take_the_diagonal_and_check_their_shape():
    pairs = np.array([[3, 0], [6, 0]])  ## W1 then W2, one pair per frame
    values = drGeometry.per_frame_distances(SOLVENT_TRAJ, pairs)
    square = md.compute_distances(SOLVENT_TRAJ, pairs, periodic=True) * 10.0
    assert square.shape == (2, 2) and not np.allclose(square[:, 0], np.diagonal(square))
    assert np.allclose(values, np.diagonal(square), atol=1e-4)
    ## blocking must not change the answer
    assert np.allclose(drGeometry.per_frame_distances(SOLVENT_TRAJ, pairs, block=1), values, atol=1e-6)
    for bad in (np.array([[3, 0]]), np.array([[3, 0, 1], [6, 0, 1]])):
        try:
            drGeometry.per_frame_distances(SOLVENT_TRAJ, bad)
        except ValueError:
            continue
        raise AssertionError(f"per_frame_distances should refuse {bad.shape}")
    try:
        drGeometry.per_frame_angles(SOLVENT_TRAJ, np.array([[3, 0, 1]]))
    except ValueError:
        return
    raise AssertionError("per_frame_angles should refuse one triplet for two frames")


def test_solvent_distances_use_the_minimum_image():
    ## the only water near C1 (at the origin) is written just inside the opposite face of the box
    waters = np.array([[[BOX - 1.5, 0.0, 0.0], [10.0, 10.0, 10.0]]])
    traj = build_system(waters)
    value = drGeometry.measure(traj, traj.topology, drGeometry.residue_labels(traj.topology, None),
                               {"name": "nearest", "type": "nearestSolvent", "a": C1})
    assert np.isclose(value[0], 1.5, atol=1e-4), f"expected the periodic image at 1.5 A, got {value[0]:.2f}"


def test_streaming_in_chunks_reproduces_the_whole_trajectory():
    ## four frames, the two designed ones plus two jittered copies
    rng = np.random.default_rng(0)
    frames = np.concatenate([WATERS, WATERS + rng.normal(scale=0.3, size=WATERS.shape)])
    traj = build_system(frames)
    topology = traj.topology
    labels = drGeometry.residue_labels(topology, None)
    specs = [{"name": "d", "type": "distance", "a": C1, "b": NE2},
             {"name": "com", "type": "comDistance", "a": {"chain": "A", "resId": 1, "atoms": "all"}, "b": NE2},
             {"name": "ang", "type": "angle", "a": O1, "b": C1, "c": NE2},
             {"name": "near", "type": "nearestSolvent", "a": C1},
             {"name": "reach", "type": "nearestSolventTo", "a": C1, "b": NE2},
             {"name": "bridge", "type": "solventBridge", "a": C1, "b": NE2, "cutoffA": 4.5, "cutoffB": 3.5},
             {"name": "bd", "type": "bridgeAngle", "a": C1, "b": C1, "c": O1}]
    oxygens = drGeometry.solvent_oxygen_indices(topology)
    for spec in specs:
        whole = drGeometry.measure(traj, topology, labels, spec, solventIndices=oxygens)
        for chunkSize in (1, 2, 3, 4):
            streamed = np.concatenate([drGeometry.measure(traj[start:start + chunkSize], topology, labels, spec,
                                                          solventIndices=oxygens)
                                       for start in range(0, traj.n_frames, chunkSize)])
            assert np.allclose(whole, streamed, atol=1e-6), f"{spec['name']} differs at chunk {chunkSize}"


def test_solvent_types_report_missing_selections():
    for spec in [{"name": "x", "type": "solventBridge", "a": C1},
                 {"name": "x", "type": "bridgeAngle", "a": C1, "b": C1}]:
        try:
            solvent_measure(spec)
        except ValueError as error:
            assert "needs" in str(error)
            continue
        raise AssertionError(f"expected a ValueError for {spec}")


def test_residue_labels_tolerate_recycled_water_numbering():
    """A solvated box reuses water residue numbers; that must not stop the solute being selectable."""
    topology = build_system(WATERS).topology
    for residue in list(topology.residues)[2:]:
        residue.resSeq = 1  ## as a PDB does once it has wrapped past 9999
    labels = drGeometry.residue_labels(topology, None)
    assert labels[("A", 1)] == 0 and labels[("B", 1)] == 1
    assert drGeometry.select_atoms(topology, labels, C1) == [0]
    ## two protein residues sharing a number is still an error worth stopping for
    clash = md.Topology()
    chain = clash.add_chain()
    for _ in range(2):
        residue = clash.add_residue("GLY", chain, resSeq=7)
        clash.add_atom("N", element.nitrogen, residue)
        clash.add_atom("CA", element.carbon, residue)
    try:
        drGeometry.residue_labels(clash, None)
    except ValueError as error:
        assert "not unique" in str(error)
        return
    raise AssertionError("two non-solvent residues with the same number should be refused")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    print("all drGeometry tests passed")
