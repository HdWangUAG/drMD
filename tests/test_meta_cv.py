"""
Checks for the collective-variable forces built by src/Surgery/drMeta.py.
Run with:  python tests/test_meta_cv.py   (or pytest tests/)
"""
import sys
from os import path as p
import numpy as np
import openmm, openmm.app as app, openmm.unit as unit

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from Surgery import drMeta

PDB = p.join(ROOT, "ExampleInputs", "Worked_Example_5_Metadynamics_of_Alanine_Dipeptide", "alanine_dipeptide.pdb")


def atom(resName, resId, atomName):
    return {"CHAIN_ID": "A", "RES_NAME": resName, "RES_ID": resId, "ATOM_NAME": atomName}


def custom(*atoms):
    return {"keyword": "custom", "customSelection": list(atoms)}


def evaluate_cvs(cvDicts):
    """Wraps each CV force in a zero-energy CustomCVForce and evaluates it on the example coordinates."""
    pdb = app.PDBFile(PDB)
    system = openmm.System()
    for a in pdb.topology.atoms():
        system.addParticle(a.element.mass)
    atomCoords = list(pdb.positions)
    cvForce = openmm.CustomCVForce("0")
    for i, cv in enumerate(cvDicts):
        cvForce.addCollectiveVariable(f"cv{i}", drMeta.gen_cv_force(cv, atomCoords, PDB))
    system.addForce(cvForce)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(pdb.positions)
    values = cvForce.getCollectiveVariableValues(context)
    return values, pdb, system


def test_com_distance_matches_distance_for_single_atoms():
    a1, a2 = atom("ACE", 1, "C"), atom("NME", 3, "N")
    values, pdb, _ = evaluate_cvs([
        {"biasVar": "distance", "selection": custom(a1, a2)},
        {"biasVar": "com_distance", "selection": custom(a1), "selection2": custom(a2)},
    ])
    print(f"DISTANCE = {values[0] * 10:.5f} A, COM_DISTANCE (one atom per group) = {values[1] * 10:.5f} A")
    assert abs(values[0] - values[1]) < 1e-9


def test_com_distance_multi_atom_groups():
    ## COM of the ACE methyl (CH3 + 3 H) to COM of the ALA side chain (CB + 3 HB)
    g1 = [atom("ACE", 1, n) for n in ["CH3", "H1", "H2", "H3"]]
    g2 = [atom("ALA", 2, n) for n in ["CB", "HB1", "HB2", "HB3"]]
    values, pdb, system = evaluate_cvs([{"biasVar": "com_distance", "selection": custom(*g1), "selection2": custom(*g2)}])
    pos = np.array(pdb.positions.value_in_unit(unit.nanometer))
    masses = np.array([system.getParticleMass(i).value_in_unit(unit.amu) for i in range(system.getNumParticles())])
    names = [(a.residue.name, int(a.residue.id), a.name) for a in pdb.topology.atoms()]
    idx1 = [names.index((d["RES_NAME"], d["RES_ID"], d["ATOM_NAME"])) for d in g1]
    idx2 = [names.index((d["RES_NAME"], d["RES_ID"], d["ATOM_NAME"])) for d in g2]
    com = lambda idx: (pos[idx] * masses[idx, None]).sum(0) / masses[idx].sum()
    expected = np.linalg.norm(com(idx1) - com(idx2))
    print(f"COM_DISTANCE = {values[0] * 10:.5f} A, numpy mass-weighted COM distance = {expected * 10:.5f} A")
    assert abs(values[0] - expected) < 1e-6


def test_torsion_angle_units():
    phi = {"biasVar": "torsion", "selection": custom(atom("ACE", 1, "C"), atom("ALA", 2, "N"), atom("ALA", 2, "CA"), atom("ALA", 2, "C"))}
    ang = {"biasVar": "angle", "selection": custom(atom("ALA", 2, "N"), atom("ALA", 2, "CA"), atom("ALA", 2, "C"))}
    values, pdb, _ = evaluate_cvs([phi, ang])
    phiDeg = values[0] * 180 / np.pi
    angDeg = values[1] * 180 / np.pi
    print(f"phi = {phiDeg:.2f} deg, N-CA-C angle = {angDeg:.2f} deg")
    assert -180 <= phiDeg <= 180 and 100 < angDeg < 125


def test_bias_variable_units_and_periodicity():
    pdb = app.PDBFile(PDB)
    atomCoords = list(pdb.positions)
    torsion = drMeta.gen_bias_variable({"biasVar": "torsion", "minValue": -180, "maxValue": 180, "biasWidth": 5.73,
                                        "selection": custom(atom("ACE", 1, "C"), atom("ALA", 2, "N"), atom("ALA", 2, "CA"), atom("ALA", 2, "C"))}, atomCoords, PDB)
    assert torsion.periodic and abs(torsion.minValue + np.pi) < 1e-9 and abs(torsion.maxValue - np.pi) < 1e-9
    dist = drMeta.gen_bias_variable({"biasVar": "com_distance", "minValue": 2.0, "maxValue": 10.0, "biasWidth": 0.5,
                                     "selection": custom(atom("ACE", 1, "C")), "selection2": custom(atom("NME", 3, "N"))}, atomCoords, PDB)
    assert not dist.periodic and abs(dist.minValue - 0.2) < 1e-12 and abs(dist.maxValue - 1.0) < 1e-12
    print("bias variable units / periodicity OK")


if __name__ == "__main__":
    test_com_distance_matches_distance_for_single_atoms()
    test_com_distance_multi_atom_groups()
    test_torsion_angle_units()
    test_bias_variable_units_and_periodicity()
    print("ALL META CV TESTS PASSED")
