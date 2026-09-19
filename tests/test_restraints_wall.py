"""
Checks for the comDistanceWall restraint (one-sided harmonic wall on a centre-of-mass distance)
built by src/Surgery/drRestraints.py, and its handling in drConfigTriage, drGaMD and drMethodsWriter.
Run with:  python tests/test_restraints_wall.py   (or pytest tests/)
"""
import sys, copy, io, contextlib
from os import path as p
import numpy as np
import openmm, openmm.app as app, openmm.unit as unit

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from Surgery import drRestraints, drGaMD
from Triage import drConfigTriage
from UtilitiesCloset import drMethodsWriter

PDB = p.join(ROOT, "ExampleInputs", "Worked_Example_5_Metadynamics_of_Alanine_Dipeptide", "alanine_dipeptide.pdb")

K = 500         ## kJ mol^-1 nm^-2
UPPER = 3.0     ## angstrom
K_NUMBER = 3


def atom(resName, resId, atomName):
    return {"CHAIN_ID": "A", "RES_NAME": resName, "RES_ID": resId, "ATOM_NAME": atomName}


def custom(*atoms):
    return {"keyword": "custom", "customSelection": list(atoms)}


## group 1: the ACE methyl (CH3 + 3 H), group 2: the NME nitrogen
GROUP1 = [atom("ACE", 1, n) for n in ["CH3", "H1", "H2", "H3"]]
GROUP2 = [atom("NME", 3, "N")]


def build_wall_system(k=K, upper=UPPER, kNumber=K_NUMBER):
    """A vacuum alanine dipeptide system whose only force is the centre-of-mass distance wall."""
    pdb = app.PDBFile(PDB)
    system = openmm.System()
    for a in pdb.topology.atoms():
        system.addParticle(a.element.mass)
    system = drRestraints.create_com_distance_wall(system, custom(*GROUP1), custom(*GROUP2), {"k": k, "upper": upper}, kNumber, PDB)
    wall = system.getForce(0)
    return pdb, system, wall


def atom_indexes(pdb, group):
    names = [(a.residue.name, int(a.residue.id), a.name) for a in pdb.topology.atoms()]
    return [names.index((d["RES_NAME"], d["RES_ID"], d["ATOM_NAME"])) for d in group]


def com_distance(pos, masses, idx1, idx2):
    com = lambda idx: (pos[idx] * masses[idx, None]).sum(0) / masses[idx].sum()
    return np.linalg.norm(com(idx1) - com(idx2))


def positions_at_com_distance(pdb, system, target_nm):
    """Translates group 2 along the COM-COM vector so that the centre-of-mass distance equals target_nm."""
    pos = np.array(pdb.positions.value_in_unit(unit.nanometer))
    masses = np.array([system.getParticleMass(i).value_in_unit(unit.amu) for i in range(system.getNumParticles())])
    idx1, idx2 = atom_indexes(pdb, GROUP1), atom_indexes(pdb, GROUP2)
    com = lambda idx: (pos[idx] * masses[idx, None]).sum(0) / masses[idx].sum()
    direction = com(idx2) - com(idx1)
    currentDistance = np.linalg.norm(direction)
    shift = direction / currentDistance * (target_nm - currentDistance)
    pos[idx2] += shift
    assert abs(com_distance(pos, masses, idx1, idx2) - target_nm) < 1e-12
    return pos


def wall_energy(system, pos):
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(pos * unit.nanometer)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)


def test_wall_energy_below_and_above_upper():
    pdb, system, wall = build_wall_system()
    upper_nm = UPPER / 10
    ## (a) exactly zero below (and at) the wall
    for d_nm in [0.5 * upper_nm, 0.9 * upper_nm, upper_nm]:
        energy = wall_energy(system, positions_at_com_distance(pdb, system, d_nm))
        assert energy == 0.0, (d_nm, energy)
    ## (b) 0.5 * k * (d - upper)^2 above the wall, k in kJ/mol/nm^2 and d in nm
    for d_nm in [1.1 * upper_nm, 2 * upper_nm, upper_nm + 0.35]:
        energy = wall_energy(system, positions_at_com_distance(pdb, system, d_nm))
        expected = 0.5 * K * (d_nm - upper_nm) ** 2
        print(f"COM distance {d_nm * 10:.3f} A: wall energy {energy:.6f} kJ/mol, expected {expected:.6f} kJ/mol")
        assert abs(energy - expected) / expected < 1e-6, (energy, expected)
    print("wall energy below / above upper OK")


def test_wall_parameters_and_groups():
    pdb, system, wall = build_wall_system()
    ## (c) k<N> is the only global parameter, upper is the only per-bond parameter
    assert isinstance(wall, openmm.CustomCentroidBondForce)
    assert wall.getNumGlobalParameters() == 1 and wall.getGlobalParameterName(0) == f"k{K_NUMBER}"
    assert abs(wall.getGlobalParameterDefaultValue(0) - K) < 1e-12
    assert wall.getNumPerBondParameters() == 1 and wall.getPerBondParameterName(0) == "upper"
    assert wall.getNumBonds() == 1 and wall.getNumGroups() == 2
    groups, bondParameters = wall.getBondParameters(0)
    assert list(groups) == [0, 1] and abs(bondParameters[0] - UPPER / 10) < 1e-12
    ## groups hold the selected atoms (mass-weighted by default: no explicit weights)
    assert list(wall.getGroupParameters(0)[0]) == atom_indexes(pdb, GROUP1)
    assert list(wall.getGroupParameters(1)[0]) == atom_indexes(pdb, GROUP2)
    assert len(wall.getGroupParameters(0)[1]) == 0
    print("wall global / per-bond parameters and groups OK")


def test_restraints_handler_wiring():
    """restraints_handler builds the wall from a restraintInfo entry and numbers k after the other restraints."""
    pdb = app.PDBFile(PDB)
    system = openmm.System()
    for a in pdb.topology.atoms():
        system.addParticle(a.element.mass)
    sim = {"restraintInfo": [
        {"restraintType": "distance", "parameters": {"k": 100, "r0": 3}, "selection": custom(atom("ACE", 1, "C"), atom("NME", 3, "N"))},
        {"restraintType": "comDistanceWall", "parameters": {"k": K, "upper": UPPER}, "selection": custom(*GROUP1), "selection2": custom(*GROUP2)},
    ]}
    system = drRestraints.restraints_handler(system, None, None, sim, None, PDB)
    forces = system.getForces()
    assert isinstance(forces[0], openmm.CustomBondForce) and isinstance(forces[1], openmm.CustomCentroidBondForce)
    assert forces[1].getGlobalParameterName(0) == "k1"
    print("restraints_handler wiring OK")


def test_clear_all_restraints_strips_k_parameters():
    ## (d) a small state XML with a wall (k1) and a position restraint (k0) among ordinary parameters
    import tempfile, os
    xml = ('<?xml version="1.0" ?>\n<State openmmVersion="8.1" time="1.0" type="State" version="1">\n'
           '\t<Parameters k0="1000" k1="500" someOtherParameter="2.5"/>\n'
           '\t<Positions>\n\t\t<Position x="0" y="0" z="0"/>\n\t</Positions>\n</State>\n')
    with tempfile.TemporaryDirectory() as tmpDir:
        xmlFile = p.join(tmpDir, "state.xml")
        with open(xmlFile, "w") as f:
            f.write(xml)
        drRestraints.clear_all_restraints(xmlFile)
        with open(xmlFile) as f:
            cleaned = f.read()
        assert not os.path.exists(xmlFile + ".tmp")
    assert 'k0=' not in cleaned and 'k1=' not in cleaned
    assert '<Parameters someOtherParameter="2.5"/>' in cleaned, cleaned
    assert '<Position x="0" y="0" z="0"/>' in cleaned
    print("clear_all_restraints strips k<N> OK")


def wall_entry(**overrides):
    entry = {"restraintType": "comDistanceWall", "parameters": {"k": 500, "upper": 14},
             "selection": {"keyword": "custom", "customSelection": [{"CHAIN_ID": "C", "RES_NAME": "S12", "RES_ID": 36, "ATOM_NAME": ["C10", "C11", "C12"]}]},
             "selection2": {"keyword": "custom", "customSelection": [{"CHAIN_ID": "A", "RES_NAME": "THR", "RES_ID": 54, "ATOM_NAME": ["CB", "OG1", "CG2"]},
                                                                      {"CHAIN_ID": "A", "RES_NAME": "ILE", "RES_ID": 63, "ATOM_NAME": ["CD1"]}]}}
    entry = copy.deepcopy(entry)
    for key, value in overrides.items():
        if value is None:
            entry.pop(key)
        else:
            entry[key] = value
    return entry


def base_config(restraintInfo):
    step = {"stepName": "01_npt", "simulationType": "NPT", "duration": "1 ns", "timestep": "2 fs", "temperature": 300,
            "logInterval": "10 ps", "restraintInfo": restraintInfo}
    return {"pathInfo": {"inputDir": ROOT, "outputDir": p.join(ROOT, "outputs")},
            "hardwareInfo": {"parallelCPU": 1, "platform": "CPU", "subprocessCpus": 1},
            "miscInfo": {"pH": 7, "firstAidMaxRetries": 1, "boxGeometry": "cubic", "writeMyMethodsSection": False,
                         "skipPdbTriage": True, "trajectorySelections": [{"selection": {"keyword": "all"}}]},
            "simulationInfo": [step]}


def validate(config):
    """Returns the validated config, or None if validation failed (drSplash calls exit(1) after printing)."""
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            return drConfigTriage.validate_config(copy.deepcopy(config))
        except SystemExit:
            return None


def test_config_triage_validation():
    ## (e) the README example is accepted ...
    assert validate(base_config([wall_entry()])) is not None
    ## ... and the wall is fine next to the other restraint types
    position = {"restraintType": "position", "parameters": {"k": 1000}, "selection": {"keyword": "protein"}}
    assert validate(base_config([position, wall_entry()])) is not None
    ## missing upper, negative k, missing selection2
    assert validate(base_config([wall_entry(parameters={"k": 500})])) is None
    assert validate(base_config([wall_entry(parameters={"k": -500, "upper": 14})])) is None
    assert validate(base_config([wall_entry(selection2=None)])) is None
    ## also: non-positive upper, non-numeric upper, bad selection2
    assert validate(base_config([wall_entry(parameters={"k": 500, "upper": 0})])) is None
    assert validate(base_config([wall_entry(parameters={"k": 500, "upper": "14"})])) is None
    assert validate(base_config([wall_entry(selection2={"keyword": "custom", "customSelection": [{"CHAIN_ID": "A"}]})])) is None
    ## the parameter checker itself
    assert drConfigTriage.check_restraint_parameters("comDistanceWall", {"k": 500, "upper": 14}) == []
    problems = drConfigTriage.check_restraint_parameters("comDistanceWall", {"k": 500})
    assert problems and "upper" in problems[0]
    ## a bad first restraint must not break reporting of a good second one
    disorders, ok = drConfigTriage.check_restraintInfo([wall_entry(parameters={"k": 500}), wall_entry()], {})
    assert not ok and disorders["restraintInfo"]["restraint_0"] and disorders["restraintInfo"]["restraint_1"] is None
    print("config triage validation OK")


def test_gamd_recognises_wall_as_restraint():
    ## (f) the wall is kept out of the GaMD boost like the other restraints
    pdb, system, wall = build_wall_system()
    assert drGaMD.is_restraint_force(wall)
    plainCentroid = openmm.CustomCentroidBondForce(2, "distance(g1,g2)")
    assert not drGaMD.is_restraint_force(plainCentroid)
    print("drGaMD.is_restraint_force OK")


def test_methods_writer_text():
    import inflect
    drMethodsWriter.inflecter = inflect.engine()
    text = drMethodsWriter.get_restraints_methods_text({"restraintInfo": [wall_entry()]})
    assert text.startswith("A one-sided harmonic wall with a force constant of 500 kJ mol<sup>-1</sup> nm<sup>-2</sup> was applied to the distance between the centres of mass of")
    assert "acting only beyond 14 Å" in text and "THR54" in text and "S1236" in text
    assert drMethodsWriter.get_force_constant_units("comDistanceWall") == "kJ mol<sup>-1</sup> nm<sup>-2</sup>"
    assert drMethodsWriter.get_restraint_target(wall_entry()) == " acting only beyond 14 Å"
    print("methods text:", text)


if __name__ == "__main__":
    test_wall_energy_below_and_above_upper()
    test_wall_parameters_and_groups()
    test_restraints_handler_wiring()
    test_clear_all_restraints_strips_k_parameters()
    test_config_triage_validation()
    test_gamd_recognises_wall_as_restraint()
    test_methods_writer_text()
    print("ALL RESTRAINT WALL TESTS PASSED")
