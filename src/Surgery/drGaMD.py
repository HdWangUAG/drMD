## BASIC PYTHON LIBRARIES
import os
from os import path as p
import json
import math

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm
import  openmm.unit  as unit

## drMD LIBRARIES
from Surgery import drSim, drMeta, drFirstAid
from ExaminationRoom import drLogger, drCheckup
from ExaminationRoom.drCVReporter import CVReporter, CV_UNITS
from UtilitiesCloset import drSelector, drFixer

## CLEAN CODE
from typing import Dict, List, Tuple, Optional
from UtilitiesCloset.drCustomClasses import FilePath, DirectoryPath

########################################################################################################
## force groups used by the GaMD integrator
## everything that is not listed below stays in group 0 (bonded, non-bonded, barostat, CMMotionRemover)
FORCE_GROUP_TOTAL: int = 0          ## boosted by the total-potential boost
FORCE_GROUP_DIHEDRAL: int = 1       ## PeriodicTorsionForce: boosted by the dihedral boost (and the total boost)
FORCE_GROUP_RESTRAINT: int = 2      ## drMD restraints: never boosted
FORCE_GROUP_CV: int = 3             ## zero-energy CustomCVForce used to monitor collective variables

## GaMD stages
STAGES: List[str] = ["cmd_stats", "gamd_equil", "gamd_prod"]
## boost "channels": P = total potential, D = dihedral potential
BOOST_CHANNELS: List[str] = ["P", "D"]

KJ_PER_KCAL: float = 4.184
########################################################################################################
########################################################################################################
@drLogger.monitor_progress_decorator()
@drFirstAid.firstAid_handler()
@drCheckup.check_up_handler()
def run_gamd(prmtop: app.AmberPrmtopFile,
              inpcrd: app.AmberInpcrdFile,
                sim: Dict,
                  saveFile: FilePath,
                    outDir: DirectoryPath,
                      platform: openmm.Platform,
                        refPdb: FilePath,
                          config: Dict) -> FilePath:
    """
    Run one stage of a Gaussian accelerated MD (GaMD) protocol [Miao, Feher & McCammon, JCTC 2015].

    GaMD adds a harmonic boost potential dV(r) = 0.5 * k * (E - V(r))^2 whenever the potential energy
    V is below a threshold E, flattening the energy landscape. The boost parameters are derived from
    the statistics (Vmax, Vmin, Vavg, sigmaV) of the potential energy, so a GaMD protocol has three
    stages, each expressed as one drMD simulation step with simulationType "GAMD":

        cmd_stats   conventional MD (no boost), potential energy statistics are collected over the whole stage
        gamd_equil  boost switched on; every updateInterval steps the boost parameters are re-derived from
                    the running Vmax / Vmin and the Vavg / sigmaV of the last updateInterval steps (the
                    AMBER "ntave" scheme: the boosted ensemble's energies replace the cmd ones as they arrive,
                    so parameters converge instead of mixing two ensembles in one accumulator)
        gamd_prod   boost parameters frozen, this is the stage to analyse

    Statistics and boost parameters are persisted in <stepName>_gamd.json in the step directory:
    each stage starts from the json of the preceding GAMD step, and a resumed step reloads its own json.
    Every logInterval, the boost potentials are written to gamd.log for reweighting (see drReweight).

    Args:
        prmtop (app.AmberPrmtopFile): The topology of the system.
        inpcrd (app.AmberInpcrdFile): The coordinates of the system.
        sim (dict): The simulation parameters.
        saveFile (str): The path to the checkpoint or XML file.
        outDir (str): The path to the output directory.
        platform (openmm.Platform): The simulation platform.
        refPdb (str): The path to the reference PDB file.
        config (dict): Dictionary containing information for all simulations in this run.

    Returns:
        str: The path to the XML file containing the final state of the simulation.
    """
    stepName: str = sim["stepName"]
    protName: str = config["proteinInfo"]["proteinName"]
    drLogger.log_info(f"Running GaMD Step: {stepName} for {protName}", True)
    ## make a simulation directory
    simDir: DirectoryPath = p.join(outDir, stepName)
    os.makedirs(simDir, exist_ok=True)

    sim = drSim.process_sim_data(sim)
    gamdInfo: Dict = sim["gamdInfo"]
    stage: str = gamdInfo["stage"]

    # Create the system (PME, HBond constraints, heavy protons, restraints and barostat)
    system: openmm.System = drSim.build_system(prmtop, inpcrd, sim, saveFile, refPdb)
    ## put the dihedral and restraint forces in their own force groups so the integrator can boost them separately
    assign_force_groups(system, gamdInfo["excludeRestraintsFromBoost"])
    ## add a zero-energy force to monitor any collective variables requested in the config
    cvForce: Optional[openmm.CustomCVForce] = add_cv_monitor_force(system, gamdInfo, prmtop, inpcrd, refPdb)

    ## create the GaMD integrator and the simulation
    integrator: openmm.CustomIntegrator = make_gamd_integrator(sim["temperature"], sim["timestep"])
    simulation, integrator = drSim.build_simulation(prmtop, system, sim, platform, integrator)
    # Load state from previous simulation (or continue from checkpoint)
    simulation: app.Simulation = drSim.load_simulation_state(simulation, saveFile)

    ## get GaMD statistics and boost parameters: from this step's json if we are resuming,
    ## otherwise from the preceding GAMD step (or fresh for a cmd_stats step)
    gamdJson: FilePath = p.join(simDir, f"{stepName}_gamd.json")
    gamdState, nStepsDone = init_gamd_state(sim, config, simDir, saveFile, gamdJson)
    ## only a plain checkpoint resume gets its remaining steps from the json (firstAid has already adjusted nSteps)
    nStepsToRun: int = sim["nSteps"]
    if is_checkpoint_resume(saveFile, simDir):
        nStepsToRun = max(sim["nSteps"] - nStepsDone, 0)
        drLogger.log_info(f"Resuming {stepName} from {stepName}_gamd.json: {nStepsDone} steps already done, {nStepsToRun} to go", True)
    apply_state_to_integrator(integrator, gamdState)
    log_boost_parameters(gamdState, stage)
    ## gamd_parameters.csv records the statistics and parameters at every update, so their convergence can be checked
    parametersCsv: FilePath = p.join(simDir, "gamd_parameters.csv")
    write_parameters_row(parametersCsv, nStepsDone, gamdState)
    updateInterval: int = gamdInfo["updateInterval"]
    if stage == "gamd_equil":
        if not is_checkpoint_resume(saveFile, simDir):
            ## Vavg / sigmaV are re-estimated from the boosted ensemble one window at a time
            reset_window_statistics(integrator)
        if updateInterval * sim["timestep"].value_in_unit(unit.picoseconds) < 20:
            drLogger.log_info(f"WARNING: updateInterval of {updateInterval} steps is a short window for estimating Vavg and sigmaV; "
                              f"GaMD boost parameters may fluctuate. 50 000 steps (100 ps at 2 fs) is typical", True, True)

    # Set up reporters
    totalSteps: int = simulation.currentStep + nStepsToRun
    reportInterval: int = sim["logInterval"]
    simulation: app.Simulation = drSim.init_reporters(simDir=simDir,
                                nSteps=totalSteps,
                                reportInterval=reportInterval,
                                simulation=simulation,
                                dcdAtomSelections= config["miscInfo"]["trajectorySelections"],
                                refPdb=refPdb)
    ## gamd.log is the input to reweighting, so it is always written. drMD resets the step counter
    ## on resume, so the reporters are offset by the steps this step has already completed
    timeOffset: float = nStepsDone * sim["timestep"].value_in_unit(unit.picoseconds)
    gamdLog: FilePath = p.join(simDir, "gamd.log")
    simulation.reporters.append(GaMDLogReporter(gamdLog, reportInterval, integrator,
                                                stepOffset=nStepsDone, timeOffset=timeOffset))
    ## cv.csv is written if collective variables were requested
    if cvForce is not None:
        cvUnits: List[str] = [CV_UNITS[cv["biasVar"].upper()] for cv in gamdInfo["cvs"]]
        cvCsv: FilePath = p.join(simDir, "cv.csv")
        simulation.reporters.append(CVReporter(cvCsv, reportInterval, cvForce, cvUnits,
                                               stepOffset=nStepsDone, timeOffset=timeOffset))

    ## the json is written on the same schedule as the checkpoint, so a resume finds statistics that
    ## match the coordinates it restarts from
    simulation.reporters.append(GaMDStateReporter(gamdJson, reportInterval, integrator, gamdState, stepOffset=nStepsDone))

    ## run the simulation in chunks of updateInterval steps; during gamd_equil the boost parameters
    ## are re-derived between chunks from the running Vmax / Vmin and the Vavg / sigmaV of the last chunk
    stepsRemaining: int = nStepsToRun
    while stepsRemaining > 0:
        chunk: int = min(updateInterval, stepsRemaining)
        simulation.step(chunk)
        stepsRemaining -= chunk
        nStepsDone += chunk
        if stage == "gamd_equil":
            gamdState = read_statistics_from_integrator(integrator, gamdState)
            gamdState = update_boost_parameters(gamdState, gamdInfo)
            apply_parameters_to_integrator(integrator, gamdState)
            reset_window_statistics(integrator)
            write_parameters_row(parametersCsv, nStepsDone, gamdState)
    gamdState = read_statistics_from_integrator(integrator, gamdState)
    gamdState["stepsCompleted"] = nStepsDone
    write_gamd_json(gamdJson, gamdState)

    ## a cmd_stats step derives the parameters the next stage will start from, from its final statistics
    ## (they are not applied to this step's integrator). gamd_prod keeps the frozen parameters it ran with.
    if stage == "cmd_stats":
        gamdState = update_boost_parameters(gamdState, gamdInfo)
        write_gamd_json(gamdJson, gamdState)
    log_boost_parameters(gamdState, f"{stage} (final)")
    if stage != "gamd_equil":
        write_parameters_row(parametersCsv, nStepsDone, gamdState)

    # save result as pdb - reset chain and residue Ids
    endPointPdb: FilePath = p.join(simDir, f"{protName}.pdb")
    drSim.write_pdb(endPointPdb, simulation)
    drFixer.reset_chains_residues(refPdb, endPointPdb, config)

    ## create a PDB file with the same atoms as the trajectory
    trajectoryPdb: FilePath = p.join(simDir, "trajectory.pdb")
    drSelector.slice_pdb_file(config["miscInfo"]["trajectorySelections"], endPointPdb, trajectoryPdb)
    # save simulation as XML
    saveXml: FilePath = p.join(simDir, f"{stepName}.xml")
    simulation.saveState(saveXml)

    return saveXml
########################################################################################################
########################################################################################################
def assign_force_groups(system: openmm.System, excludeRestraintsFromBoost: bool) -> None:
    """
    Moves the PeriodicTorsionForce into FORCE_GROUP_DIHEDRAL and (optionally) drMD's restraint
    forces into FORCE_GROUP_RESTRAINT so that the GaMD integrator can read their energies and
    forces separately. Restraint forces are recognised as the Custom*Force objects created by
    drRestraints (they all carry a global parameter starting with "k").

    Args:
        system (openmm.System): the system to modify in place
        excludeRestraintsFromBoost (bool): if True, restraints are moved out of the boosted groups
    """
    for force in system.getForces():
        if isinstance(force, openmm.PeriodicTorsionForce):
            force.setForceGroup(FORCE_GROUP_DIHEDRAL)
        elif excludeRestraintsFromBoost and is_restraint_force(force):
            force.setForceGroup(FORCE_GROUP_RESTRAINT)
########################################################################################################
def is_restraint_force(force: openmm.Force) -> bool:
    """
    Decides whether a force was added by drRestraints.
    drRestraints only makes CustomExternalForce, CustomBondForce, CustomAngleForce and CustomTorsionForce
    objects, each with a global parameter named k<N>.
    """
    if not isinstance(force, (openmm.CustomExternalForce, openmm.CustomBondForce,
                              openmm.CustomAngleForce, openmm.CustomTorsionForce)):
        return False
    for i in range(force.getNumGlobalParameters()):
        if force.getGlobalParameterName(i).startswith("k"):
            return True
    return False
########################################################################################################
def add_cv_monitor_force(system: openmm.System,
                          gamdInfo: Dict,
                            prmtop: app.AmberPrmtopFile,
                              inpcrd: app.AmberInpcrdFile,
                                refPdb: FilePath) -> Optional[openmm.CustomCVForce]:
    """
    Adds a zero-energy CustomCVForce wrapping each collective variable listed in gamdInfo["cvs"].
    The force contributes nothing to the dynamics; it exists so that CVReporter can evaluate the
    variables with the same machinery (drMeta.gen_cv_force) as metadynamics.

    Returns:
        the CustomCVForce, or None if no collective variables were requested
    """
    cvs: List[Dict] = gamdInfo.get("cvs", [])
    if len(cvs) == 0:
        return None
    atomCoords: list = drMeta.get_atom_coords_for_metadynamics(prmtop, inpcrd)
    cvForce: openmm.CustomCVForce = openmm.CustomCVForce("0")
    for cvIndex, cv in enumerate(cvs):
        cvForce.addCollectiveVariable(f"cv{cvIndex}", drMeta.gen_cv_force(cv, atomCoords, refPdb))
    cvForce.setForceGroup(FORCE_GROUP_CV)
    system.addForce(cvForce)
    return cvForce
########################################################################################################
def make_gamd_integrator(temperature: unit.Quantity,
                          timestep: unit.Quantity,
                            friction: unit.Quantity = 1/unit.picosecond) -> openmm.CustomIntegrator:
    """
    Builds a Langevin-middle integrator (the same scheme as openmm.LangevinMiddleIntegrator) in which
    the force used for the velocity update is the GaMD-boosted force

        F' = F_boostable * (1 - kP * (EP - V))  -  kD * (ED - VD) * F_dihedral  +  F_restraint

    with each boost term active only while its potential is below its threshold. V is the energy of
    force groups 0 and 1, VD that of group 1, and group 2 (restraints) is never boosted.

    The integrator also accumulates the running maximum, minimum, mean and variance (Welford's
    algorithm) of V and VD every step, in global variables. Global variables are evaluated in double
    precision on every platform, so the accumulators are safe in single-precision CUDA runs.

    Global variables (all energies in kJ/mol):
        kP, EP, kD, ED             boost parameters, k in 1/(kJ/mol). k = 0 switches a boost off.
        VT, VD                     boostable potential energy, dihedral potential energy (current step)
        dVP, dVD                   boost potentials added this step
        nStats                     number of steps accumulated in the mean / variance statistics
        seeded                     0 until the first step has initialised Vmax / Vmin, then 1
        VmaxP, VminP, VavgP, M2P   running statistics of VT (M2 = sum of squared deviations)
        VmaxD, VminD, VavgD, M2D   running statistics of VD

    Args:
        temperature (unit.Quantity or int): simulation temperature (K)
        timestep (unit.Quantity): integration timestep
        friction (unit.Quantity): Langevin friction coefficient

    Returns:
        openmm.CustomIntegrator: the integrator
    """
    if not unit.is_quantity(temperature):
        temperature = temperature * unit.kelvin
    kT: float = (unit.MOLAR_GAS_CONSTANT_R * temperature).value_in_unit(unit.kilojoules_per_mole)
    dt: float = timestep.value_in_unit(unit.picoseconds)
    gamma: float = friction.value_in_unit(unit.picosecond**-1)

    integrator: openmm.CustomIntegrator = openmm.CustomIntegrator(timestep)
    ## Langevin middle constants
    integrator.addGlobalVariable("a", math.exp(-gamma * dt))
    integrator.addGlobalVariable("b", math.sqrt(1 - math.exp(-2 * gamma * dt)))
    integrator.addGlobalVariable("kT", kT)
    integrator.addPerDofVariable("x1", 0)
    integrator.addPerDofVariable("fBoost", 0)
    ## boost parameters
    for name in ["kP", "EP", "kD", "ED"]:
        integrator.addGlobalVariable(name, 0)
    ## per-step energies, boosts and force weights
    for name in ["VT", "VD", "dVP", "dVD", "fwP", "fwD", "deltaP", "deltaD"]:
        integrator.addGlobalVariable(name, 0)
    ## running statistics; "seeded" is 0 until the first energy has initialised Vmax / Vmin
    integrator.addGlobalVariable("nStats", 0)
    integrator.addGlobalVariable("seeded", 0)
    for channel in BOOST_CHANNELS:
        for name in ["Vmax", "Vmin", "Vavg", "M2"]:
            integrator.addGlobalVariable(f"{name}{channel}", 0)

    integrator.addUpdateContextState()
    ## energies of the boostable (groups 0 + 1) and dihedral (group 1) potentials
    ## NB. a single computation step may only refer to one force group, hence the two-step sums
    integrator.addComputeGlobal("VD", f"energy{FORCE_GROUP_DIHEDRAL}")
    integrator.addComputeGlobal("VT", f"energy{FORCE_GROUP_TOTAL}")
    integrator.addComputeGlobal("VT", "VT + VD")
    ## running statistics (Welford). On the very first step, max and min are seeded with the energy
    integrator.addComputeGlobal("nStats", "nStats + 1")
    for channel, energyName in zip(BOOST_CHANNELS, ["VT", "VD"]):
        integrator.addComputeGlobal(f"Vmax{channel}", f"select(seeded, max(Vmax{channel}, {energyName}), {energyName})")
        integrator.addComputeGlobal(f"Vmin{channel}", f"select(seeded, min(Vmin{channel}, {energyName}), {energyName})")
        integrator.addComputeGlobal(f"delta{channel}", f"{energyName} - Vavg{channel}")
        integrator.addComputeGlobal(f"Vavg{channel}", f"Vavg{channel} + delta{channel}/nStats")
        integrator.addComputeGlobal(f"M2{channel}", f"M2{channel} + delta{channel}*({energyName} - Vavg{channel})")
    integrator.addComputeGlobal("seeded", "1")
    ## boost potentials and the resulting force weights
    integrator.addComputeGlobal("dVP", "0.5*kP*(EP-VT)^2*step(EP-VT)")
    integrator.addComputeGlobal("dVD", "0.5*kD*(ED-VD)^2*step(ED-VD)")
    integrator.addComputeGlobal("fwP", "1 - kP*(EP-VT)*step(EP-VT)")
    integrator.addComputeGlobal("fwD", "-kD*(ED-VD)*step(ED-VD)")
    ## assemble the boosted force one force group at a time:
    ## F' = (f0 + f1)*fwP + f1*fwD + f2
    integrator.addComputePerDof("fBoost", f"f{FORCE_GROUP_TOTAL}*fwP")
    integrator.addComputePerDof("fBoost", f"fBoost + f{FORCE_GROUP_DIHEDRAL}*(fwP + fwD)")
    integrator.addComputePerDof("fBoost", f"fBoost + f{FORCE_GROUP_RESTRAINT}")
    ## Langevin middle scheme with the boosted force
    integrator.addComputePerDof("v", "v + dt*fBoost/m")
    integrator.addConstrainVelocities()
    integrator.addComputePerDof("x", "x + 0.5*dt*v")
    integrator.addComputePerDof("v", "a*v + b*sqrt(kT/m)*gaussian")
    integrator.addComputePerDof("x", "x + 0.5*dt*v")
    integrator.addComputePerDof("x1", "x")
    integrator.addConstrainPositions()
    integrator.addComputePerDof("v", "v + (x-x1)/dt")

    return integrator
########################################################################################################
def new_gamd_state(sim: Dict) -> Dict:
    """
    An empty GaMD state: no statistics, boosts switched off. All energies are stored in kcal/mol
    (the AMBER convention) and force constants k in 1/(kcal/mol); k0 is dimensionless.
    """
    gamdInfo: Dict = sim["gamdInfo"]
    return {"stepName": sim["stepName"],
            "stage": gamdInfo["stage"],
            "boostType": gamdInfo["boostType"],
            "thresholdMode": gamdInfo["thresholdMode"],
            "units": "kcal/mol",
            "extremesSeeded": False,
            "stepsCompleted": 0,
            "totalSteps": sim["nSteps"],
            "statistics": {channel: {"count": 0, "Vmax": 0.0, "Vmin": 0.0, "Vavg": 0.0, "sigmaV": 0.0}
                           for channel in BOOST_CHANNELS},
            "parameters": {channel: {"k0": 0.0, "k": 0.0, "E": 0.0, "thresholdModeUsed": "off"}
                           for channel in BOOST_CHANNELS}}
########################################################################################################
def is_checkpoint_resume(saveFile: Optional[FilePath], simDir: DirectoryPath) -> bool:
    """True if saveFile is the checkpoint file that drOperator found inside this step's own directory."""
    if saveFile is None:
        return False
    return p.splitext(saveFile)[1] == ".chk" and p.dirname(p.abspath(saveFile)) == p.abspath(simDir)
########################################################################################################
def find_previous_gamd_json(sim: Dict, config: Dict, outDir: DirectoryPath) -> Optional[FilePath]:
    """
    Finds the <stepName>_gamd.json written by the closest preceding GAMD step in the config.

    Returns:
        the path to the json, or None if this is the first GAMD step
    """
    simulations: List[Dict] = config["simulationInfo"]
    stepNames: List[str] = [simulation["stepName"] for simulation in simulations]
    thisIndex: int = stepNames.index(sim["stepName"])
    for previousSim in reversed(simulations[:thisIndex]):
        if previousSim["simulationType"].upper() == "GAMD":
            previousStepName: str = previousSim["stepName"]
            return p.join(outDir, previousStepName, f"{previousStepName}_gamd.json")
    return None
########################################################################################################
def init_gamd_state(sim: Dict,
                     config: Dict,
                       simDir: DirectoryPath,
                         saveFile: Optional[FilePath],
                           gamdJson: FilePath) -> Tuple[Dict, int]:
    """
    Decides where this step's statistics and boost parameters come from:
      1. resuming (checkpoint or firstAid restart inside this step): this step's own json
      2. cmd_stats: a fresh state
      3. gamd_equil / gamd_prod: the json of the preceding GAMD step, with the parameters
         (re)derived from its statistics for this step's boostType / thresholdMode

    Returns:
        gamdState (dict), nStepsDone (int): steps of this step already completed (only non-zero when resuming)
    """
    stage: str = sim["gamdInfo"]["stage"]
    outDir: DirectoryPath = p.dirname(simDir)
    resuming: bool = saveFile is not None and p.isfile(gamdJson) and \
        p.abspath(saveFile).startswith(p.abspath(simDir) + os.sep)
    if resuming:
        gamdState: Dict = read_gamd_json(gamdJson)
        return gamdState, gamdState["stepsCompleted"]

    if stage == "cmd_stats":
        return new_gamd_state(sim), 0

    previousJson: Optional[FilePath] = find_previous_gamd_json(sim, config, outDir)
    if previousJson is None or not p.isfile(previousJson):
        raise FileNotFoundError(f"GaMD stage {stage} in step {sim['stepName']} needs a preceding GAMD step "
                                f"(cmd_stats or gamd_equil) but no <stepName>_gamd.json was found: {previousJson}")
    previousState: Dict = read_gamd_json(previousJson)
    gamdState: Dict = new_gamd_state(sim)
    gamdState["statistics"] = previousState["statistics"]
    gamdState["extremesSeeded"] = previousState.get("extremesSeeded", True)
    gamdState["derivedFrom"] = previousJson
    if stage == "gamd_prod" and previousState["stage"] != "cmd_stats":
        ## production freezes the parameters the preceding equilibration (or production replicate) used
        gamdState["parameters"] = previousState["parameters"]
    else:
        gamdState = update_boost_parameters(gamdState, sim["gamdInfo"])
    return gamdState, 0
########################################################################################################
def compute_boost_parameters(Vmax: float, Vmin: float, Vavg: float, sigmaV: float,
                              sigma0: float, thresholdMode: str) -> Tuple[float, float, float, str]:
    """
    GaMD boost parameters from potential energy statistics (Miao, Feher & McCammon, JCTC 2015, eqs 9-13).
    All energies in the same unit (kcal/mol here).

    lower bound (E = Vmax):
        k0 = min(1, (sigma0/sigmaV) * (Vmax - Vmin)/(Vmax - Vavg))
    upper bound (E = Vmin + (Vmax - Vmin)/k0):
        k0 = (1 - sigma0/sigmaV) * (Vmax - Vmin)/(Vavg - Vmin), valid only when 0 < k0 <= 1,
        otherwise falls back to the lower bound.

    Returns:
        k0 (dimensionless), k (1/energy), E (energy), thresholdModeUsed ("lower" | "upper" | "off")
    """
    energyRange: float = Vmax - Vmin
    if energyRange <= 0 or sigmaV <= 0 or Vmax - Vavg <= 0:
        return 0.0, 0.0, Vmax, "off"

    thresholdModeUsed: str = "lower"
    if thresholdMode.lower() == "upper":
        if Vavg - Vmin > 0:
            k0: float = (1 - sigma0 / sigmaV) * energyRange / (Vavg - Vmin)
            if 0 < k0 <= 1:
                E: float = Vmin + energyRange / k0
                return k0, k0 / energyRange, E, "upper"
        drLogger.log_info("GaMD: upper-bound threshold gives k0 outside (0, 1], falling back to lower bound (E = Vmax)", True, True)

    k0: float = min(1.0, (sigma0 / sigmaV) * energyRange / (Vmax - Vavg))
    E: float = Vmax
    return k0, k0 / energyRange, E, thresholdModeUsed
########################################################################################################
def update_boost_parameters(gamdState: Dict, gamdInfo: Dict) -> Dict:
    """
    (Re)derives the boost parameters for both channels from the statistics held in gamdState,
    switching off whichever channel the boostType does not use.
    """
    boostType: str = gamdInfo["boostType"].lower()
    activeChannels: List[str] = {"total": ["P"], "dihedral": ["D"], "dual": ["P", "D"]}[boostType]
    sigma0: Dict[str, float] = {"P": gamdInfo["sigma0P"], "D": gamdInfo["sigma0D"]}
    for channel in BOOST_CHANNELS:
        stats: Dict = gamdState["statistics"][channel]
        if channel in activeChannels and stats["count"] > 0:
            k0, k, E, modeUsed = compute_boost_parameters(stats["Vmax"], stats["Vmin"], stats["Vavg"], stats["sigmaV"],
                                                          sigma0[channel], gamdInfo["thresholdMode"])
        else:
            k0, k, E, modeUsed = 0.0, 0.0, stats["Vmax"], "off"
        gamdState["parameters"][channel] = {"k0": k0, "k": k, "E": E, "thresholdModeUsed": modeUsed}
    return gamdState
########################################################################################################
def apply_state_to_integrator(integrator: openmm.CustomIntegrator, gamdState: Dict) -> None:
    """
    Loads boost parameters and running statistics from gamdState (kcal/mol) into the integrator's
    global variables (kJ/mol).
    """
    apply_parameters_to_integrator(integrator, gamdState)
    for channel in BOOST_CHANNELS:
        stats: Dict = gamdState["statistics"][channel]
        integrator.setGlobalVariableByName(f"Vmax{channel}", stats["Vmax"] * KJ_PER_KCAL)
        integrator.setGlobalVariableByName(f"Vmin{channel}", stats["Vmin"] * KJ_PER_KCAL)
        integrator.setGlobalVariableByName(f"Vavg{channel}", stats["Vavg"] * KJ_PER_KCAL)
        integrator.setGlobalVariableByName(f"M2{channel}", stats["count"] * (stats["sigmaV"] * KJ_PER_KCAL) ** 2)
    ## both channels share the same step count
    integrator.setGlobalVariableByName("nStats", gamdState["statistics"]["P"]["count"])
    integrator.setGlobalVariableByName("seeded", 1 if gamdState.get("extremesSeeded", False) else 0)
########################################################################################################
def reset_window_statistics(integrator: openmm.CustomIntegrator) -> None:
    """
    Restarts the mean / variance accumulators (Vavg, M2, nStats) so that the next read gives the statistics
    of the steps run from now on. Vmax and Vmin are kept: they are running extremes over the whole protocol
    (the "seeded" flag stays set, so they are not re-initialised).
    """
    integrator.setGlobalVariableByName("nStats", 0)
    for channel in BOOST_CHANNELS:
        integrator.setGlobalVariableByName(f"Vavg{channel}", 0.0)
        integrator.setGlobalVariableByName(f"M2{channel}", 0.0)
########################################################################################################
def apply_parameters_to_integrator(integrator: openmm.CustomIntegrator, gamdState: Dict) -> None:
    """Loads only the boost parameters (kcal/mol) into the integrator (kJ/mol), leaving the running statistics untouched."""
    for channel in BOOST_CHANNELS:
        params: Dict = gamdState["parameters"][channel]
        integrator.setGlobalVariableByName(f"k{channel}", params["k"] / KJ_PER_KCAL)
        integrator.setGlobalVariableByName(f"E{channel}", params["E"] * KJ_PER_KCAL)
########################################################################################################
def read_statistics_from_integrator(integrator: openmm.CustomIntegrator, gamdState: Dict) -> Dict:
    """
    Reads the running statistics accumulated by the integrator (kJ/mol) back into gamdState (kcal/mol).
    """
    count: int = int(round(integrator.getGlobalVariableByName("nStats")))
    gamdState["extremesSeeded"] = integrator.getGlobalVariableByName("seeded") > 0.5
    for channel in BOOST_CHANNELS:
        M2: float = integrator.getGlobalVariableByName(f"M2{channel}")
        sigmaV: float = math.sqrt(max(M2, 0.0) / count) if count > 0 else 0.0
        gamdState["statistics"][channel] = {"count": count,
                                            "Vmax": integrator.getGlobalVariableByName(f"Vmax{channel}") / KJ_PER_KCAL,
                                            "Vmin": integrator.getGlobalVariableByName(f"Vmin{channel}") / KJ_PER_KCAL,
                                            "Vavg": integrator.getGlobalVariableByName(f"Vavg{channel}") / KJ_PER_KCAL,
                                            "sigmaV": sigmaV / KJ_PER_KCAL}
    return gamdState
########################################################################################################
def write_gamd_json(gamdJson: FilePath, gamdState: Dict) -> None:
    """Writes the GaMD state atomically (write to a temp file, then rename)."""
    tempJson: FilePath = gamdJson + ".tmp"
    with open(tempJson, "w") as f:
        json.dump(gamdState, f, indent=2)
    os.replace(tempJson, gamdJson)
########################################################################################################
def read_gamd_json(gamdJson: FilePath) -> Dict:
    with open(gamdJson, "r") as f:
        return json.load(f)
########################################################################################################
def write_parameters_row(parametersCsv: FilePath, step: int, gamdState: Dict) -> None:
    """
    Appends one row of statistics and boost parameters (kcal/mol) to gamd_parameters.csv.
    Columns: step, then for each channel (P, D): count, Vmax, Vmin, Vavg, sigmaV, k0, k, E, threshold.
    """
    writeHeader: bool = not (p.isfile(parametersCsv) and os.path.getsize(parametersCsv) > 0)
    with open(parametersCsv, "a") as f:
        if writeHeader:
            columns: List[str] = ["step"]
            for channel in BOOST_CHANNELS:
                columns += [f"{name}_{channel}" for name in ["count", "Vmax", "Vmin", "Vavg", "sigmaV", "k0", "k", "E", "threshold"]]
            f.write(",".join(columns) + "\n")
        row: List[str] = [str(step)]
        for channel in BOOST_CHANNELS:
            stats: Dict = gamdState["statistics"][channel]
            params: Dict = gamdState["parameters"][channel]
            row += [str(stats["count"])] + [f"{stats[name]:.4f}" for name in ["Vmax", "Vmin", "Vavg", "sigmaV"]]
            row += [f"{params['k0']:.6f}", f"{params['k']:.8f}", f"{params['E']:.4f}", params["thresholdModeUsed"]]
        f.write(",".join(row) + "\n")
########################################################################################################
def log_boost_parameters(gamdState: Dict, label: str) -> None:
    for channel, name in zip(BOOST_CHANNELS, ["total", "dihedral"]):
        stats: Dict = gamdState["statistics"][channel]
        params: Dict = gamdState["parameters"][channel]
        drLogger.log_info(f"GaMD {label} | {name:8s} | Vmax {stats['Vmax']:.2f} Vmin {stats['Vmin']:.2f} "
                          f"Vavg {stats['Vavg']:.2f} sigmaV {stats['sigmaV']:.2f} kcal/mol (n = {stats['count']}) | "
                          f"k0 {params['k0']:.4f} E {params['E']:.2f} threshold {params['thresholdModeUsed']}")
########################################################################################################
class GaMDStateReporter(object):
    """
    Writes <stepName>_gamd.json at the same interval as drMD's checkpoint reporter, holding the running
    statistics read from the integrator, the boost parameters currently applied and the number of steps
    completed, so that a resumed step can carry on exactly where the checkpoint left off.
    """
    def __init__(self, gamdJson: FilePath, reportInterval: int, integrator: openmm.CustomIntegrator,
                 gamdState: Dict, stepOffset: int = 0) -> None:
        self._gamdJson: FilePath = gamdJson
        self._reportInterval: int = reportInterval
        self._integrator: openmm.CustomIntegrator = integrator
        self._gamdState: Dict = gamdState      ## shared with run_gamd, which updates the parameters in place
        self._stepOffset: int = stepOffset

    def describeNextReport(self, simulation: app.Simulation) -> tuple:
        steps: int = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation: app.Simulation, state: openmm.State) -> None:
        read_statistics_from_integrator(self._integrator, self._gamdState)
        self._gamdState["stepsCompleted"] = simulation.currentStep + self._stepOffset
        write_gamd_json(self._gamdJson, self._gamdState)
########################################################################################################
class GaMDLogReporter(object):
    """
    Writes gamd.log: one row per logInterval with the boostable potential energy, the dihedral
    potential energy and the two boost potentials, all in kcal/mol. The energies are evaluated
    for the reported frame (the same frame the trajectory reporter writes), so rows align 1:1
    with trajectory frames and with cv.csv. This file is the only input reweighting needs.

    Columns: step, time_ps, V_total, V_dihedral, dV_P, dV_D
    """
    def __init__(self, file: FilePath, reportInterval: int, integrator: openmm.CustomIntegrator,
                 append: bool = True, stepOffset: int = 0, timeOffset: float = 0.0) -> None:
        self._reportInterval: int = reportInterval
        self._integrator: openmm.CustomIntegrator = integrator
        ## drMD resets the step counter on resume; the offsets keep step / time monotonic across segments
        self._stepOffset: int = stepOffset
        self._timeOffset: float = timeOffset
        writeHeader: bool = not (append and p.isfile(file) and os.path.getsize(file) > 0)
        self._out = open(file, "a" if append else "w")
        if writeHeader:
            self._out.write("step,time_ps,V_total,V_dihedral,dV_P,dV_D\n")
            self._out.flush()

    def describeNextReport(self, simulation: app.Simulation) -> tuple:
        steps: int = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation: app.Simulation, state: openmm.State) -> None:
        context: openmm.Context = simulation.context
        VT: float = context.getState(getEnergy=True, groups={FORCE_GROUP_TOTAL, FORCE_GROUP_DIHEDRAL}).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        VD: float = context.getState(getEnergy=True, groups={FORCE_GROUP_DIHEDRAL}).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        kP: float = self._integrator.getGlobalVariableByName("kP")
        EP: float = self._integrator.getGlobalVariableByName("EP")
        kD: float = self._integrator.getGlobalVariableByName("kD")
        ED: float = self._integrator.getGlobalVariableByName("ED")
        dVP: float = 0.5 * kP * (EP - VT) ** 2 if VT < EP else 0.0
        dVD: float = 0.5 * kD * (ED - VD) ** 2 if VD < ED else 0.0
        timePs: float = state.getTime().value_in_unit(unit.picoseconds) + self._timeOffset
        self._out.write(f"{simulation.currentStep + self._stepOffset},{timePs:.3f},{VT / KJ_PER_KCAL:.5f},{VD / KJ_PER_KCAL:.5f},"
                        f"{dVP / KJ_PER_KCAL:.5f},{dVD / KJ_PER_KCAL:.5f}\n")
        self._out.flush()

    def __del__(self) -> None:
        try:
            self._out.close()
        except Exception:
            pass
########################################################################################################
