"""
Numerical checks for the GaMD integrator in src/Surgery/drGaMD.py.
Run with:  python tests/test_gamd_integrator.py   (or pytest tests/)
Uses the alanine dipeptide from Worked Example 5, parameterised in vacuum with tleap-free
OpenMM force field files so the test needs no AmberTools.
"""
import os, sys, math
from os import path as p
import numpy as np
import openmm, openmm.app as app, openmm.unit as unit

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from Surgery import drGaMD

PDB = p.join(ROOT, "ExampleInputs", "Worked_Example_5_Metadynamics_of_Alanine_Dipeptide", "alanine_dipeptide.pdb")
KJ = unit.kilojoules_per_mole


def make_vacuum_system(restraint: bool):
    pdb = app.PDBFile(PDB)
    ff = app.ForceField("amber14-all.xml")
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, constraints=None, removeCMMotion=False)
    if restraint:
        ## the same shape of force drRestraints makes: a CustomExternalForce with global parameter k0
        rest = openmm.CustomExternalForce("k0*periodicdistance(x, y, z, x0, y0, z0)^2")
        rest.addGlobalParameter("k0", 500.0)
        for name in ["x0", "y0", "z0"]:
            rest.addPerParticleParameter(name)
        for i, pos in enumerate(modeller.positions):
            if i % 3 == 0:
                rest.addParticle(i, [pos.x + 0.05, pos.y, pos.z])  # slightly displaced so it exerts a force
        system.addForce(rest)
    return system, modeller.topology, modeller.positions


def group_forces_and_energies(context):
    out = {}
    for g in [drGaMD.FORCE_GROUP_TOTAL, drGaMD.FORCE_GROUP_DIHEDRAL, drGaMD.FORCE_GROUP_RESTRAINT]:
        state = context.getState(getForces=True, getEnergy=True, groups={g})
        out[g] = (state.getForces(asNumpy=True).value_in_unit(unit.kilojoules_per_mole / unit.nanometer),
                  state.getPotentialEnergy().value_in_unit(KJ))
    return out


def run_one_step(kP, EP_offset, kD, ED_offset, restraint, platformName="Reference"):
    """One frictionless integrator step; compare velocity update to the analytical boosted force."""
    system, topology, positions = make_vacuum_system(restraint)
    drGaMD.assign_force_groups(system, excludeRestraintsFromBoost=True)
    dt = 0.0005 * unit.picoseconds
    integ = drGaMD.make_gamd_integrator(300 * unit.kelvin, dt, friction=0 / unit.picosecond)
    platform = openmm.Platform.getPlatformByName(platformName)
    context = openmm.Context(system, integ, platform)
    context.setPositions(positions)
    rng = np.random.default_rng(0)
    v0 = rng.normal(0, 0.5, size=(system.getNumParticles(), 3))
    context.setVelocities(v0 * unit.nanometers / unit.picosecond)

    groups = group_forces_and_energies(context)
    f0, e0 = groups[0]; f1, e1 = groups[1]; f2, e2 = groups[2]
    VT = e0 + e1; VD = e1
    EP = VT + EP_offset; ED = VD + ED_offset
    for name, val in [("kP", kP), ("EP", EP), ("kD", kD), ("ED", ED)]:
        integ.setGlobalVariableByName(name, val)

    integ.step(1)

    ## analytical boosted force
    fwP = 1 - kP * (EP - VT) if VT < EP else 1.0
    fwD = -kD * (ED - VD) if VD < ED else 0.0
    fBoost = (f0 + f1) * fwP + f1 * fwD + f2
    masses = np.array([system.getParticleMass(i).value_in_unit(unit.amu) for i in range(system.getNumParticles())])
    v1_expected = v0 + dt.value_in_unit(unit.picoseconds) * fBoost / masses[:, None]
    v1 = context.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(unit.nanometers / unit.picosecond)
    err = np.abs(v1 - v1_expected).max()
    scale = np.abs(v1_expected - v0).max()

    ## integrator-side bookkeeping
    dVP_expected = 0.5 * kP * (EP - VT) ** 2 if VT < EP else 0.0
    dVD_expected = 0.5 * kD * (ED - VD) ** 2 if VD < ED else 0.0
    got = {n: integ.getGlobalVariableByName(n) for n in ["VT", "VD", "dVP", "dVD", "nStats", "VmaxP", "VminP", "VavgP", "M2P"]}
    assert abs(got["VT"] - VT) < 1e-6 * max(1, abs(VT)), (got["VT"], VT)
    assert abs(got["VD"] - VD) < 1e-6 * max(1, abs(VD)), (got["VD"], VD)
    assert abs(got["dVP"] - dVP_expected) < 1e-8 * max(1, abs(dVP_expected)), (got["dVP"], dVP_expected)
    assert abs(got["dVD"] - dVD_expected) < 1e-8 * max(1, abs(dVD_expected)), (got["dVD"], dVD_expected)
    assert got["nStats"] == 1 and abs(got["VmaxP"] - VT) < 1e-6 and abs(got["VminP"] - VT) < 1e-6 and abs(got["VavgP"] - VT) < 1e-6 and got["M2P"] == 0
    return err, scale, e2


def test_force_expression():
    cases = [  # kP, EP offset, kD, ED offset, restraint
        (0.0, 0.0, 0.0, 0.0, False),          # no boost: plain leapfrog step
        (0.01, 50.0, 0.0, 0.0, False),        # total boost only (fwP = 0.5)
        (0.0, 0.0, 0.5, 20.0, False),         # dihedral boost only
        (0.01, 50.0, 0.5, 20.0, False),       # dual boost
        (0.01, 50.0, 0.5, 20.0, True),        # dual boost with an (unboosted) restraint force present
        (0.01, -50.0, 0.5, -20.0, True),      # thresholds below the energy: boost must switch off
    ]
    for case in cases:
        err, scale, e2 = run_one_step(*case)
        rel = err / scale
        print(f"case kP={case[0]} EP-V={case[1]} kD={case[2]} ED-VD={case[3]} restraint={case[4]}: "
              f"max |v1 - expected| = {err:.2e} (relative {rel:.2e}), restraint energy {e2:.2f} kJ/mol")
        assert rel < 1e-6, f"boosted force mismatch: relative error {rel}"
        if case[4]:
            assert e2 > 0, "restraint force should be present in group 2"


def test_welford_statistics():
    """Statistics accumulated by the integrator must match numpy over a short trajectory."""
    system, topology, positions = make_vacuum_system(False)
    drGaMD.assign_force_groups(system, True)
    integ = drGaMD.make_gamd_integrator(300 * unit.kelvin, 0.001 * unit.picoseconds)
    context = openmm.Context(system, integ, openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    context.setVelocitiesToTemperature(300 * unit.kelvin, 1)
    energies = []
    for _ in range(200):
        integ.step(1)
        energies.append(integ.getGlobalVariableByName("VT"))
    energies = np.array(energies)
    n = integ.getGlobalVariableByName("nStats")
    assert n == 200
    assert abs(integ.getGlobalVariableByName("VmaxP") - energies.max()) < 1e-6
    assert abs(integ.getGlobalVariableByName("VminP") - energies.min()) < 1e-6
    assert abs(integ.getGlobalVariableByName("VavgP") - energies.mean()) < 1e-6
    sigma = math.sqrt(integ.getGlobalVariableByName("M2P") / n)
    assert abs(sigma - energies.std()) < 1e-6 * max(1, energies.std())
    print(f"Welford statistics over 200 steps agree with numpy: mean {energies.mean():.3f}, std {energies.std():.3f} kJ/mol")


def test_boost_parameters():
    """Boost parameter formulas (Miao 2015) and the upper-bound fallback."""
    k0, k, E, mode = drGaMD.compute_boost_parameters(Vmax=-1000.0, Vmin=-1200.0, Vavg=-1100.0, sigmaV=10.0, sigma0=6.0, thresholdMode="lower")
    assert mode == "lower" and E == -1000.0
    assert abs(k0 - min(1, (6 / 10) * 200 / 100)) < 1e-12 and abs(k - k0 / 200) < 1e-12
    ## upper bound valid: k0 = (1 - 6/10) * 200/100 = 0.8
    k0, k, E, mode = drGaMD.compute_boost_parameters(-1000.0, -1200.0, -1100.0, 10.0, 6.0, "upper")
    assert mode == "upper" and abs(k0 - 0.8) < 1e-12 and abs(E - (-1200 + 200 / 0.8)) < 1e-9
    ## upper bound invalid (sigma0 > sigmaV -> k0 < 0): falls back to lower
    k0, k, E, mode = drGaMD.compute_boost_parameters(-1000.0, -1200.0, -1100.0, 4.0, 6.0, "upper")
    assert mode == "lower" and E == -1000.0 and k0 == 1.0
    ## degenerate statistics switch the boost off
    k0, k, E, mode = drGaMD.compute_boost_parameters(-1000.0, -1000.0, -1000.0, 0.0, 6.0, "lower")
    assert mode == "off" and k0 == 0.0 and k == 0.0
    print("boost parameter formulas OK")


if __name__ == "__main__":
    test_boost_parameters()
    test_welford_statistics()
    test_force_expression()
    print("ALL GaMD INTEGRATOR TESTS PASSED")
