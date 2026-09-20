"""
Checks for the optional flat bottom on distance and torsion restraints built by
src/Surgery/drRestraints.py, and its handling in drConfigTriage, drGaMD and drMethodsWriter.
A flat-bottomed restraint applies no force while the measured value is within halfWidth of
the target, and a harmonic penalty in the excess beyond it.
Run with:  python tests/test_restraints_flat_bottom.py   (or pytest tests/)
"""
import sys, math
from os import path as p
import numpy as np
import openmm, openmm.app as app, openmm.unit as unit

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from Surgery import drRestraints, drGaMD
from Triage import drConfigTriage
from UtilitiesCloset import drMethodsWriter

PDB = p.join(ROOT, "ExampleInputs", "Worked_Example_5_Metadynamics_of_Alanine_Dipeptide", "alanine_dipeptide.pdb")

K_DIST = 1000.0    ## kJ mol^-1 nm^-2
R0 = 3.4           ## angstrom      (the build spec's Asp281 Od - thioester C1 target)
HW_DIST = 0.4      ## angstrom      (its half-width)
K_TORS = 500.0     ## kJ mol^-1 rad^-2
PHI0 = -70.0       ## degrees       (the family His285 chi1)
HW_TORS = 20.0     ## degrees       (its half-width)


def atom(resName, resId, atomName):
    return {"CHAIN_ID": "A", "RES_NAME": resName, "RES_ID": resId, "ATOM_NAME": atomName}


def custom(*atoms):
    return {"keyword": "custom", "customSelection": list(atoms)}


PAIR = [atom("ACE", 1, "CH3"), atom("NME", 3, "N")]
QUARTET = [atom("ACE", 1, "C"), atom("ALA", 2, "N"), atom("ALA", 2, "CA"), atom("ALA", 2, "C")]


def empty_system(pdb):
    system = openmm.System()
    for a in pdb.topology.atoms():
        system.addParticle(a.element.mass)
    return system


def atom_indexes(pdb, group):
    names = [(a.residue.name, int(a.residue.id), a.name) for a in pdb.topology.atoms()]
    return [names.index((d["RES_NAME"], d["RES_ID"], d["ATOM_NAME"])) for d in group]


def context_for(system):
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    return openmm.Context(system, integrator, openmm.Platform.getPlatformByName("Reference"))


def energy_and_forces(context, positions_nm):
    context.setPositions(positions_nm * unit.nanometer)
    state = context.getState(getEnergy=True, getForces=True)
    return (state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
            np.array(state.getForces().value_in_unit(unit.kilojoule_per_mole / unit.nanometer)))


def place_pair(pdb, idx, separation_nm):
    """All atoms at the origin except the restrained pair, which is put on the x axis."""
    pos = np.zeros((pdb.topology.getNumAtoms(), 3))
    pos[idx[0]] = [0.0, 0.0, 0.0]
    pos[idx[1]] = [separation_nm, 0.0, 0.0]
    return pos


def dihedral(p0, p1, p2, p3):
    """The i-j-k-l dihedral in degrees, IUPAC sign convention."""
    b0, axis, b2 = p0 - p1, p2 - p1, p3 - p2
    axis = axis / np.linalg.norm(axis)
    v = b0 - np.dot(b0, axis) * axis
    w = b2 - np.dot(b2, axis) * axis
    return math.degrees(math.atan2(np.dot(np.cross(axis, v), w), np.dot(v, w)))


def place_dihedral(pdb, idx, phi_degrees):
    """Four atoms in a standard dihedral frame; the fourth is rotated to the requested angle."""
    phi = math.radians(phi_degrees)
    pos = np.zeros((pdb.topology.getNumAtoms(), 3))
    pos[idx[0]] = [0.0, 0.1, 0.0]                                            ## i, off the j-k axis
    pos[idx[1]] = [0.0, 0.0, 0.0]                                            ## j
    pos[idx[2]] = [0.15, 0.0, 0.0]                                           ## k, the axis runs along x
    pos[idx[3]] = [0.15, 0.1 * math.cos(phi), 0.1 * math.sin(phi)]           ## l, rotated about j-k
    placed = dihedral(*(pos[i] for i in idx))
    assert abs((placed - phi_degrees + 180) % 360 - 180) < 1e-6, f"placed {placed}°, wanted {phi_degrees}°"
    return pos


def expected_flat_bottom(value, target, halfWidth, k, periodic=False):
    deviation = value - target
    if periodic:   ## torsions: wrap into (-pi, pi], so 190 degrees away is really 170 the other way
        deviation = (deviation + math.pi) % (2 * math.pi) - math.pi
    return 0.5 * k * max(0.0, abs(deviation) - halfWidth) ** 2


def test_distance_flat_bottom():
    pdb = app.PDBFile(PDB)
    system = empty_system(pdb)
    system = drRestraints.create_distance_restraint(
        system, custom(*PAIR), {"k": K_DIST, "r0": R0, "halfWidth": HW_DIST}, 0, PDB)
    idx = atom_indexes(pdb, PAIR)
    context = context_for(system)

    for separation_angstrom in [R0, R0 - HW_DIST, R0 + HW_DIST, R0 - 0.2, R0 + 0.3,
                                R0 - 1.0, R0 + 1.0, R0 + 2.5, 0.5]:
        energy, forces = energy_and_forces(context, place_pair(pdb, idx, separation_angstrom / 10.0))
        ## the expression works in nm, so convert the angstrom deviation before squaring
        expected = expected_flat_bottom(separation_angstrom / 10.0, R0 / 10.0, HW_DIST / 10.0, K_DIST)
        assert abs(energy - expected) < 1e-6, f"r = {separation_angstrom} Å: {energy} != {expected}"
        if abs(separation_angstrom - R0) < HW_DIST:
            assert np.abs(forces).max() < 1e-6, f"force inside the flat bottom at r = {separation_angstrom} Å"
        elif abs(separation_angstrom - R0) > HW_DIST + 0.05:
            assert np.abs(forces[idx]).max() > 1e-3, f"no force outside the flat bottom at r = {separation_angstrom} Å"
    print("  distance flat bottom: zero inside ±%.1f Å of %.1f Å, harmonic in the excess outside" % (HW_DIST, R0))


def test_distance_without_half_bottom_is_unchanged():
    """No halfWidth (or zero) must reproduce the original plain harmonic restraint."""
    pdb = app.PDBFile(PDB)
    idx = atom_indexes(pdb, PAIR)
    for parameters in [{"k": K_DIST, "r0": R0}, {"k": K_DIST, "r0": R0, "halfWidth": 0}]:
        system = drRestraints.create_distance_restraint(empty_system(pdb), custom(*PAIR), parameters, 0, PDB)
        context = context_for(system)
        for separation_angstrom in [R0, R0 - 1.0, R0 + 1.0]:
            energy, _ = energy_and_forces(context, place_pair(pdb, idx, separation_angstrom / 10.0))
            expected = 0.5 * K_DIST * ((separation_angstrom - R0) / 10.0) ** 2
            assert abs(energy - expected) < 1e-6, f"{parameters}: {energy} != {expected}"
    print("  distance without a flat bottom: unchanged plain harmonic")


def test_torsion_flat_bottom():
    pdb = app.PDBFile(PDB)
    system = empty_system(pdb)
    system = drRestraints.create_torsion_restraint(
        system, custom(*QUARTET), {"k": K_TORS, "phi0": PHI0, "halfWidth": HW_TORS}, 0, PDB)
    idx = atom_indexes(pdb, QUARTET)
    context = context_for(system)

    for phi in [PHI0, PHI0 - HW_TORS, PHI0 + HW_TORS, PHI0 - 10, PHI0 + 15,
                PHI0 - 45, PHI0 + 45, 0.0, 120.0]:
        energy, forces = energy_and_forces(context, place_dihedral(pdb, idx, phi))
        expected = expected_flat_bottom(math.radians(phi), math.radians(PHI0), math.radians(HW_TORS), K_TORS, periodic=True)
        assert abs(energy - expected) < 1e-5, f"phi = {phi}°: {energy} != {expected}"
        if abs(phi - PHI0) < HW_TORS - 1:
            assert np.abs(forces).max() < 1e-5, f"force inside the flat bottom at phi = {phi}°"
    print("  torsion flat bottom: zero inside ±%.0f° of %.0f°, harmonic in the excess outside" % (HW_TORS, PHI0))


def test_torsion_flat_bottom_is_periodic():
    """A target near +180° must see -170° as 10° away, not 350° away."""
    pdb = app.PDBFile(PDB)
    system = drRestraints.create_torsion_restraint(
        empty_system(pdb), custom(*QUARTET), {"k": K_TORS, "phi0": 175.0, "halfWidth": 5.0}, 0, PDB)
    idx = atom_indexes(pdb, QUARTET)
    context = context_for(system)

    energy, _ = energy_and_forces(context, place_dihedral(pdb, idx, -178.0))   ## 7° away
    expected = 0.5 * K_TORS * math.radians(2.0) ** 2
    assert abs(energy - expected) < 1e-5, f"wrapped deviation: {energy} != {expected}"
    energy, _ = energy_and_forces(context, place_dihedral(pdb, idx, -179.0))   ## 6° away
    assert abs(energy - 0.5 * K_TORS * math.radians(1.0) ** 2) < 1e-5
    energy, _ = energy_and_forces(context, place_dihedral(pdb, idx, 178.0))    ## inside the well
    assert energy < 1e-6, f"wrapped flat bottom: {energy}"
    print("  torsion flat bottom: deviation wrapped into (-180, 180]")


def test_gamd_sees_the_flat_bottom_forces_as_restraints():
    """drGaMD must keep flat-bottomed restraints out of the boosted force groups."""
    pdb = app.PDBFile(PDB)
    system = empty_system(pdb)
    system = drRestraints.create_distance_restraint(
        system, custom(*PAIR), {"k": K_DIST, "r0": R0, "halfWidth": HW_DIST}, 0, PDB)
    system = drRestraints.create_torsion_restraint(
        system, custom(*QUARTET), {"k": K_TORS, "phi0": PHI0, "halfWidth": HW_TORS}, 1, PDB)
    ## the plain harmonic forms must be recognised too - they hold k per bond / per torsion,
    ## not as a global parameter, so they were being left in the boosted group
    system = drRestraints.create_distance_restraint(system, custom(*PAIR), {"k": K_DIST, "r0": R0}, 2, PDB)
    system = drRestraints.create_torsion_restraint(system, custom(*QUARTET), {"k": K_TORS, "phi0": PHI0}, 3, PDB)
    drGaMD.assign_force_groups(system, excludeRestraintsFromBoost=True)
    for force in system.getForces():
        assert force.getForceGroup() == drGaMD.FORCE_GROUP_RESTRAINT, (force.__class__.__name__, force.getForceGroup())
    assert drGaMD.is_restraint_force(openmm.HarmonicBondForce()) is False
    print("  drGaMD: flat-bottomed and plain restraints put in the unboosted restraint force group")


def test_config_triage_accepts_and_rejects_half_width():
    ok = drConfigTriage.check_restraint_parameters("distance", {"k": K_DIST, "r0": R0, "halfWidth": HW_DIST})
    assert ok == [], ok
    ok = drConfigTriage.check_restraint_parameters("torsion", {"k": K_TORS, "phi0": PHI0, "halfWidth": HW_TORS})
    assert ok == [], ok
    ## no half-width at all is still valid
    assert drConfigTriage.check_restraint_parameters("distance", {"k": K_DIST, "r0": R0}) == []
    ## and the ways it can be wrong
    for restraintType, parameters, expect in [
            ("distance", {"k": K_DIST, "r0": R0, "halfWidth": -1}, "negative"),
            ("distance", {"k": K_DIST, "r0": R0, "halfWidth": "wide"}, "number"),
            ("torsion", {"k": K_TORS, "phi0": PHI0, "halfWidth": 200}, "less than 180"),
            ("position", {"k": 10, "halfWidth": 1}, "only supported"),
            ("angle", {"k": 10, "theta0": 90, "halfWidth": 1}, "only supported")]:
        problems = drConfigTriage.check_restraint_parameters(restraintType, parameters)
        assert any(expect in problem for problem in problems), (restraintType, parameters, problems)
    print("  drConfigTriage: halfWidth checked on distance and torsion, refused elsewhere")


def test_methods_writer_mentions_the_flat_bottom():
    sim = {"restraintInfo": [
        {"restraintType": "distance", "parameters": {"k": K_DIST, "r0": R0, "halfWidth": HW_DIST},
         "selection": custom(*PAIR)},
        {"restraintType": "torsion", "parameters": {"k": K_TORS, "phi0": PHI0, "halfWidth": HW_TORS},
         "selection": custom(*QUARTET)},
        {"restraintType": "distance", "parameters": {"k": K_DIST, "r0": R0}, "selection": custom(*PAIR)}]}
    drMethodsWriter.inflecter = drMethodsWriter.inflect.engine()   ## normally set up by the writer's entry point
    text = drMethodsWriter.get_restraints_methods_text(sim)
    assert text.count("flat-bottomed") == 2, text
    assert f"within {HW_DIST} Å of the target" in text, text
    assert f"within {HW_TORS} degrees of the target" in text, text
    print("  drMethodsWriter: flat bottom described once per flat-bottomed restraint")


if __name__ == "__main__":
    for test in [test_distance_flat_bottom, test_distance_without_half_bottom_is_unchanged,
                 test_torsion_flat_bottom, test_torsion_flat_bottom_is_periodic,
                 test_gamd_sees_the_flat_bottom_forces_as_restraints,
                 test_config_triage_accepts_and_rejects_half_width,
                 test_methods_writer_mentions_the_flat_bottom]:
        print(test.__name__)
        test()
    print("\nall flat-bottom restraint tests passed")
