## BASIC PYTHON LIBRARIES
import os
from os import path as p

## OPENMM LIBRARIES
import openmm.app as app
import openmm as openmm
import openmm.unit as unit

## CLEAN CODE
from typing import List, Optional
from UtilitiesCloset.drCustomClasses import FilePath

#####################################################################################
## unit conversions from OpenMM internal units (nm, radians) to the units used in drMD configs
CV_UNIT_FACTORS = {"angstrom": 10.0,
                   "degrees": 180.0 / 3.141592653589793,
                   "none": 1.0}
## which unit each collective variable type is reported in
CV_UNITS = {"TORSION": "degrees",
            "ANGLE": "degrees",
            "DISTANCE": "angstrom",
            "COM_DISTANCE": "angstrom",
            "RMSD": "angstrom"}
#####################################################################################
class CVReporter(object):
    """
    An OpenMM reporter that writes the values of collective variables (and, optionally,
    the energy of a bias force acting on them) to a CSV file at a fixed interval.

    Follows the same reporter protocol as app.StateDataReporter, so it can be appended to
    simulation.reporters alongside the reporters made in drSim.init_reporters.

    Columns: step, time_ps, cv_0 ... cv_n, bias_energy_kJ
    Collective variables are written in drMD's config units (angstrom for distances / RMSD,
    degrees for angles / torsions). The bias energy is in kJ/mol.
    """
    def __init__(self,
                  file: FilePath,
                    reportInterval: int,
                      cvForce: openmm.CustomCVForce,
                        cvUnits: List[str],
                          biasForceGroup: Optional[int] = None,
                            append: bool = True,
                              stepOffset: int = 0,
                                timeOffset: float = 0.0) -> None:
        """
        Args:
            file (FilePath): path to the CSV file to write
            reportInterval (int): the interval (in time steps) at which to write a row
            cvForce (openmm.CustomCVForce): the force whose collective variables are reported.
                For metadynamics this is the bias force created by openmm.app.metadynamics.Metadynamics,
                for GaMD it is a zero-energy CustomCVForce used purely to evaluate the variables.
            cvUnits (List[str]): one of "angstrom", "degrees" or "none" for each collective variable
            biasForceGroup (int, optional): force group of the bias force. If given, the potential
                energy of that group is written in the bias_energy_kJ column, otherwise 0 is written.
            append (bool): open the file in append mode so a resumed simulation continues the same file
            stepOffset (int): added to the step column. drMD resets the step counter when a simulation is
                resumed, so pass the number of steps completed before the resume to keep steps monotonic.
            timeOffset (float): the same for the time_ps column (ps)
        """
        self._reportInterval: int = reportInterval
        self._stepOffset: int = stepOffset
        self._timeOffset: float = timeOffset
        self._cvForce: openmm.CustomCVForce = cvForce
        self._factors: List[float] = [CV_UNIT_FACTORS[cvUnit] for cvUnit in cvUnits]
        self._biasForceGroup: Optional[int] = biasForceGroup
        nCvs: int = cvForce.getNumCollectiveVariables()
        ## only write a header if we are starting a fresh file
        writeHeader: bool = not (append and p.isfile(file) and os.path.getsize(file) > 0)
        self._out = open(file, "a" if append else "w")
        if writeHeader:
            cvNames: List[str] = [f"cv_{i}" for i in range(nCvs)]
            self._out.write(",".join(["step", "time_ps"] + cvNames + ["bias_energy_kJ"]) + "\n")
            self._out.flush()

    def describeNextReport(self, simulation: app.Simulation) -> tuple:
        steps: int = self._reportInterval - simulation.currentStep % self._reportInterval
        ## (steps, positions, velocities, forces, energies, enforcePeriodicBox)
        return (steps, False, False, False, False, None)

    def report(self, simulation: app.Simulation, state: openmm.State) -> None:
        cvValues: list = self._cvForce.getCollectiveVariableValues(simulation.context)
        cvValues = [value * factor for value, factor in zip(cvValues, self._factors)]
        if self._biasForceGroup is None:
            biasEnergy: float = 0.0
        else:
            biasState: openmm.State = simulation.context.getState(getEnergy=True, groups={self._biasForceGroup})
            biasEnergy: float = biasState.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        timePs: float = state.getTime().value_in_unit(unit.picoseconds) + self._timeOffset
        row: List[str] = [str(simulation.currentStep + self._stepOffset), f"{timePs:.3f}"] + [f"{value:.5f}" for value in cvValues] + [f"{biasEnergy:.5f}"]
        self._out.write(",".join(row) + "\n")
        self._out.flush()

    def __del__(self) -> None:
        try:
            self._out.close()
        except Exception:
            pass
#####################################################################################
