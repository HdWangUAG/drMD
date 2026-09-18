## BASIC PYTHON LIBRARIES
import os
from os import path as p
import numpy as np
import pandas as pd

## OPENMM LIBRARIES
import openmm.app as app
from openmm.app import metadynamics
import openmm as openmm
import  openmm.unit  as unit

## drMD LIBRARIES
from Surgery import drSim, drFirstAid
from ExaminationRoom import drLogger, drCheckup
from UtilitiesCloset import drSelector, drFixer

########################################################################################################
########################################################################################################
@drLogger.monitor_progress_decorator()
@drFirstAid.firstAid_handler()
@drCheckup.check_up_handler()
def run_metadynamics(prmtop: app.Topology,
                      inpcrd: any,
                        sim: dict,
                          saveFile: str,
                            outDir: str,
                              platform: openmm.Platform,
                                refPdb: str,
                                    config:dict) -> str:

    """
    Run a simulation at constant pressure (NpT) step with biases.

    Args:
        prmtop (str): The path to the topology file.
        inpcrd (str): The path to the coordinates file.
        sim (dict): The simulation parameters.
        saveFile (str): The path to the checkpoint or XML file.
        simDir (str): The path to the simulation directory.
        platform (openmm.Platform): The simulation platform.
        pdbFile (str): The path to the PDB file.

    Returns:
        str: The path to the XML file containing the final state of the simulation.

    This function runs a metadynamics simulation at constant pressure (NpT) step with biases.
    It initializes the system, handles any constraints, adds a constant pressure force,
    sets up the integrator and system, loads the state from a checkpoint or XML file,
    sets up reporters, runs the simulation, saves the final geometry as a PDB file,
    resets the chain and residue Ids, runs drCheckup, performs trajectory clustering
    if specified in the simulation parameters, and saves the simulation state as an
    XML file.
    """
    stepName = sim["stepName"]
    drLogger.log_info(f"Running MetaDynamics Step: {stepName}",True)
    ## make a simulation directory
    simDir: str = p.join(outDir, stepName)
    os.makedirs(simDir, exist_ok=True)

    sim = drSim.process_sim_data(sim)
    # Create the system (PME, HBond constraints, heavy protons, restraints and barostat)
    system: openmm.System = drSim.build_system(prmtop, inpcrd, sim, saveFile, refPdb)

    # Read metaDynamicsInfo from sim config
    metaDynamicsInfo: dict = sim["metaDynamicsInfo"]

    # Read biases from sim config and create bias variables
    biases: list = metaDynamicsInfo["biases"]
    biasVariables: list = []

    atomCoords: list = get_atom_coords_for_metadynamics(prmtop, inpcrd)
    for bias in biases:
        # Create bias variable based on the type of bias
        biasVariable: metadynamics.BiasVariable = gen_bias_variable(bias, atomCoords, refPdb)
        biasVariables.append(biasVariable)

    # Create metadynamics object and add bias variables as forces to the system
    meta: metadynamics.Metadynamics = metadynamics.Metadynamics(system=system,
                                     variables=biasVariables,
                                     temperature=sim["temperature"],
                                     biasFactor=metaDynamicsInfo["biasFactor"],
                                     height=metaDynamicsInfo["height"],
                                     frequency=50,
                                     saveFrequency=50,
                                     biasDir=simDir)
    

    # Set up integrator and create new simulation
    simulation, integrator = drSim.build_simulation(prmtop, system, sim, platform)
    # Load state from previous simulation (or continue from checkpoint)
    simulation: app.Simulation = drSim.load_simulation_state(simulation, saveFile)
    # Set up reporters
    totalSteps: int = simulation.currentStep + sim["nSteps"]
    reportInterval: int = sim["logInterval"]
    simulation: app.Simulation = drSim.init_reporters(simDir=simDir,
                                nSteps=totalSteps,
                                reportInterval=reportInterval,
                                simulation=simulation,
                                dcdAtomSelections= config["miscInfo"]["trajectorySelections"],
                                refPdb=refPdb)
    # Run metadynamics simulation
    meta.step(simulation, sim["nSteps"])
 
    # find name to call outFiles
    protName: str = p.basename(p.dirname(simDir))
    # save result as pdb - reset chain and residue Ids
    endPointPdb: str = p.join(simDir, f"{protName}.pdb")

    drSim.write_pdb(endPointPdb, simulation)
    # with open(endPointPdb, 'w') as output:
    #     app.pdbfile.PDBFile.writeFile(simulation.topology, state.getPositions(), output)
    ## reset the chain and residue ids to original
    drFixer.reset_chains_residues(refPdb, endPointPdb, config)

    ## create a PDB file with the same atoms as the trajectory
    trajectoryPdb = p.join(simDir, "trajectory.pdb")
    drSelector.slice_pdb_file(config["miscInfo"]["trajectorySelections"], endPointPdb, trajectoryPdb)
    # save simulation as XML
    saveXml: str = p.join(simDir, f"{stepName}.xml")
    simulation.saveState(saveXml)

    ## get free energy and write to csv
    metadynamicsFreeEnergy = meta.getFreeEnergy()
    freeEnergyCsv = p.join(simDir, "freeEnergy.csv")
    np.savetxt(freeEnergyCsv, metadynamicsFreeEnergy, delimiter=",", header="Bias Variable Free Energy")

    # Return checkpoint file for continuing simulation
    return saveXml
########################################################################################################
########################################################################################################
def get_atom_coords_for_metadynamics(prmtop: app.Topology, inpcrd: any) -> list:
    """
    Get the coordinates of all atoms in the system.

    Parameters
    ----------
    prmtop : file
        The topology file.
    inpcrd : file
        The coordinates file.

    Returns
    -------
    atomCoords : list
        The coordinates of all atoms in the system.
    """
    # Get the topology and positions from the input files
    topology: app.Topology = prmtop.topology
    positions: np.ndarray = inpcrd.positions

    # Initialize an empty list to store atom coordinates
    atomCoords: list = []

    # Loop over all atoms in the topology
    for i, atom in enumerate(topology.atoms()):
        # Append the coordinates of the current atom to the atomCoords list
        atomCoords.append(positions[i])

    # Return the list of atom coordinates
    return atomCoords
########################################################################################################
def gen_cv_force(bias: dict, atomCoords: list, refPdb: str) -> openmm.Force:
    """
    Create the openmm.Force whose energy defines a collective variable.
    Used both for metadynamics bias variables and for monitoring collective variables
    during other simulation types (e.g. GaMD).

    Parameters
    ----------
    bias : dict
        A bias / collective-variable dictionary from the config. Must contain "biasVar" and "selection";
        COM_DISTANCE additionally needs "selection2".
    atomCoords : list
        The coordinates of all atoms in the system (used as the RMSD reference).
    refPdb : str
        The reference PDB file used to resolve selections.

    Returns
    -------
    cvForce : openmm.Force
        One of RMSDForce, CustomTorsionForce, CustomBondForce, CustomAngleForce or CustomCentroidBondForce.
    """
    biasVar: str = bias["biasVar"].upper()
    atomIndexes: list = drSelector.get_atom_indexes(bias["selection"], refPdb)

    if biasVar == "RMSD":
        ## RMSD to the starting coordinates of the selected atoms
        cvForce: openmm.RMSDForce = openmm.RMSDForce(atomCoords, atomIndexes)

    elif biasVar == "TORSION":
        ## dihedral angle between four atoms
        cvForce: openmm.CustomTorsionForce = openmm.CustomTorsionForce("theta")
        cvForce.addTorsion(atomIndexes[0],
                            atomIndexes[1],
                              atomIndexes[2],
                                atomIndexes[3])

    elif biasVar == "DISTANCE":
        ## distance between two atoms
        cvForce: openmm.CustomBondForce = openmm.CustomBondForce("r")
        cvForce.addBond(atomIndexes[0],
                         atomIndexes[1])

    elif biasVar == "COM_DISTANCE":
        ## distance between the (mass-weighted) centres of mass of two groups of atoms
        ## group 1 comes from "selection", group 2 from "selection2"
        atomIndexes2: list = drSelector.get_atom_indexes(bias["selection2"], refPdb)
        cvForce: openmm.CustomCentroidBondForce = openmm.CustomCentroidBondForce(2, "distance(g1,g2)")
        cvForce.addGroup(atomIndexes)
        cvForce.addGroup(atomIndexes2)
        cvForce.addBond([0, 1])

    elif biasVar == "ANGLE":
        ## angle between three atoms
        cvForce: openmm.CustomAngleForce = openmm.CustomAngleForce("theta")
        cvForce.addAngle(atomIndexes[0],
                          atomIndexes[1],
                            atomIndexes[2])
    else:
        raise ValueError(f"Unknown biasVar {bias['biasVar']}, must be one of RMSD, TORSION, DISTANCE, COM_DISTANCE, ANGLE")

    return cvForce
########################################################################################################
def gen_bias_variable(bias: dict, atomCoords: list, refPdb: str) -> metadynamics.BiasVariable:
    """
    Generate a metadynamics bias variable from a bias dictionary.
    Distances and RMSDs are given in angstrom, angles and torsions in degrees.
    Torsions are periodic, everything else is not.

    Parameters
    ----------
    bias : dict
        The bias dictionary containing biasVar, minValue, maxValue, biasWidth and selection(s).
    atomCoords : list
        The coordinates of all atoms in the system.
    refPdb : str
        The reference PDB file used to resolve selections.

    Returns
    -------
    biasVariable : metadynamics.BiasVariable
        The bias variable.
    """
    biasVar: str = bias["biasVar"].upper()
    cvForce: openmm.Force = gen_cv_force(bias, atomCoords, refPdb)

    if biasVar in ["RMSD", "DISTANCE", "COM_DISTANCE"]:
        cvUnit: unit.Unit = unit.angstrom
    else:
        cvUnit: unit.Unit = unit.degrees

    biasVariable: metadynamics.BiasVariable = metadynamics.BiasVariable(force = cvForce,
                                                    minValue = bias["minValue"] * cvUnit,
                                                    maxValue = bias["maxValue"] * cvUnit,
                                                    biasWidth = bias["biasWidth"] * cvUnit,
                                                    periodic = biasVar == "TORSION")
    return biasVariable
########################################################################################################
