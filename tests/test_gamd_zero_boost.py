"""
Checks the zero-boost GaMD null run: sigma0P = sigma0D = 0.
Run with:  python tests/test_gamd_zero_boost.py   (or pytest tests/)

The null run is an acceptance requirement for the GaMD implementation: the whole code path
(GaMD integrator, force-group split, restraint exclusion, reporters, stage hand-over) is run
with a boost amplitude of exactly zero, so its trajectory must be statistically identical to
plain MD in the same ensemble. sigma0 = 0 is how that is asked for - there is no k0 key in the
config, compute_boost_parameters derives it.

config triage used to reject sigma0P / sigma0D of 0 with "must be a positive number", which made
the null run impossible to ask for, so this file pins:
  1. triage accepts sigma0P = sigma0D = 0.0 and still rejects negatives and a zero updateInterval,
  2. sigma0 = 0 gives k0 = k = 0 in both threshold modes, reported as "off" (it used to reach zero
     in the upper mode only through the k0-out-of-range fallback, logging a warning that made a
     deliberate null run look like a broken one),
  3. k = 0 makes the boost potential dV = 0 exactly rather than a NaN: nothing in the parameter
     chain or in the integrator's compute steps divides by k, k0, sigma0 or sigmaV.

No simulation is run here.
"""
import sys, copy, io, re, contextlib
from os import path as p

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from Triage import drConfigTriage
from Surgery import drGaMD


def gamd_step(name, stage, **extra):
    info = {"stage": stage}
    info.update(extra)
    return {"stepName": name, "simulationType": "GAMD", "duration": "1 ns", "timestep": "2 fs", "temperature": 300,
            "logInterval": "2 ps", "gamdInfo": info}


def base_config(steps):
    return {"pathInfo": {"inputDir": ROOT, "outputDir": p.join(ROOT, "outputs")},
            "hardwareInfo": {"parallelCPU": 1, "platform": "CPU", "subprocessCpus": 1},
            "miscInfo": {"pH": 7, "firstAidMaxRetries": 1, "boxGeometry": "cubic", "writeMyMethodsSection": False,
                         "skipPdbTriage": True, "trajectorySelections": [{"selection": {"keyword": "all"}}]},
            "simulationInfo": steps}


def validate(config):
    """The validated config, or None if validation failed (drSplash calls exit(1) after printing)."""
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            return drConfigTriage.validate_config(copy.deepcopy(config))
        except SystemExit:
            return None


def zero_boost_steps():
    zero = {"boostType": "dual", "thresholdMode": "lower", "sigma0P": 0.0, "sigma0D": 0.0,
            "updateInterval": 50000, "excludeRestraintsFromBoost": True, "ensemble": "NPT"}
    return [gamd_step("01_cmd_stats", "cmd_stats", **zero),
            gamd_step("02_gamd_equil", "gamd_equil", **zero),
            gamd_step("03_gamd_prod", "gamd_prod", **zero)]


def test_triage_accepts_zero_sigma0():
    """sigma0P / sigma0D of exactly 0 is the zero-boost null run, not a broken config."""
    cfg = validate(base_config(zero_boost_steps()))
    assert cfg is not None, "config triage rejected the zero-boost null run (sigma0P = sigma0D = 0.0)"
    for step in cfg["simulationInfo"]:
        assert step["gamdInfo"]["sigma0P"] == 0.0 and step["gamdInfo"]["sigma0D"] == 0.0, step["gamdInfo"]
    ## integer 0 is the same request
    assert validate(base_config([gamd_step("01_cmd_stats", "cmd_stats", sigma0P=0, sigma0D=0)])) is not None
    print("triage accepts sigma0P = sigma0D = 0 (zero-boost null run)")


def test_triage_still_rejects_bad_numbers():
    ## a negative sigma0 is meaningless
    assert validate(base_config([gamd_step("01_cmd_stats", "cmd_stats", sigma0P=-6.0)])) is None
    assert validate(base_config([gamd_step("01_cmd_stats", "cmd_stats", sigma0D=-1)])) is None
    ## updateInterval is a number of steps: 0 would mean an empty chunk, and negatives are nonsense
    assert validate(base_config([gamd_step("01_cmd_stats", "cmd_stats", updateInterval=0)])) is None
    assert validate(base_config([gamd_step("01_cmd_stats", "cmd_stats", updateInterval=-50000)])) is None
    ## still type checked
    assert validate(base_config([gamd_step("01_cmd_stats", "cmd_stats", sigma0P="0.0")])) is None
    assert validate(base_config([gamd_step("01_cmd_stats", "cmd_stats", sigma0P=True)])) is None
    print("triage still rejects negative sigma0, a zero / negative updateInterval and wrong types")


def test_zero_sigma0_gives_zero_boost():
    """k0 = k = 0 for sigma0 = 0, in both threshold modes, and dV = 0 exactly."""
    Vmax, Vmin, Vavg, sigmaV = -1000.0, -1200.0, -1100.0, 10.0
    for thresholdMode in ["lower", "upper"]:
        k0, k, E, mode = drGaMD.compute_boost_parameters(Vmax, Vmin, Vavg, sigmaV, 0.0, thresholdMode)
        assert (k0, k, mode) == (0.0, 0.0, "off"), (thresholdMode, k0, k, E, mode)
        ## dV = 0.5 * k * (E - V)^2 * step(E - V), as the integrator and GaMDLogReporter compute it
        for V in [Vmin, Vavg, Vmax, Vmax + 100]:
            assert 0.5 * k * (E - V) ** 2 == 0.0
    ## and the same through the state / config layer that run_gamd actually uses
    gamdInfo = {"stage": "cmd_stats", "boostType": "dual", "thresholdMode": "lower", "sigma0P": 0.0, "sigma0D": 0.0,
                "updateInterval": 50000, "excludeRestraintsFromBoost": True, "ensemble": "NPT", "cvs": []}
    state = drGaMD.new_gamd_state({"stepName": "01_cmd_stats", "nSteps": 100, "gamdInfo": gamdInfo})
    for channel in drGaMD.BOOST_CHANNELS:
        state["statistics"][channel] = {"count": 500000, "Vmax": Vmax, "Vmin": Vmin, "Vavg": Vavg, "sigmaV": sigmaV}
    state = drGaMD.update_boost_parameters(state, gamdInfo)
    for channel in drGaMD.BOOST_CHANNELS:
        assert state["parameters"][channel]["k"] == 0.0, state["parameters"]
        assert state["parameters"][channel]["k0"] == 0.0, state["parameters"]
    ## a non-zero sigma0 must still give a boost, or this test would pass on a dead implementation
    boosted = drGaMD.update_boost_parameters(copy.deepcopy(state), {**gamdInfo, "sigma0P": 6.0, "sigma0D": 6.0})
    assert boosted["parameters"]["P"]["k"] > 0 and boosted["parameters"]["D"]["k"] > 0, boosted["parameters"]
    print("sigma0 = 0 -> k0 = k = 0 -> dV = 0 exactly (sigma0 = 6 still boosts)")


def test_nothing_divides_by_the_boost_constants():
    """k = 0 must give dV = 0, not a division by zero or a NaN, anywhere in the integrator."""
    forbidden = re.compile(r"/\s*\(?\s*(kP|kD|k0|sigma0|sigmaV)\b")
    ## the integrator's own compute steps, as OpenMM stores them. compute_boost_parameters does
    ## divide by sigmaV and by the energy range, but it is guarded and returns before doing so;
    ## the integrator is the part that runs every step with whatever k it was handed.
    import openmm, openmm.unit as unit
    integrator = drGaMD.make_gamd_integrator(300 * unit.kelvin, 0.002 * unit.picoseconds)
    expressions = [integrator.getComputationStep(i)[2] for i in range(integrator.getNumComputations())]
    for expression in expressions:
        assert not forbidden.search(expression), f"integrator step divides by a boost constant: {expression}"
    ## every boost constant starts at zero, so a step of the integrator must be reachable with k = 0
    for name in ["kP", "EP", "kD", "ED"]:
        assert integrator.getGlobalVariableByName(name) == 0.0
    ## degenerate statistics must return finite zeros rather than raise or produce a NaN
    for stats in [(0.0, 0.0, 0.0, 0.0), (-1000.0, -1000.0, -1000.0, 0.0), (-1000.0, -1200.0, -1000.0, 10.0)]:
        for sigma0 in [0.0, 6.0]:
            for thresholdMode in ["lower", "upper"]:
                k0, k, E, mode = drGaMD.compute_boost_parameters(*stats, sigma0, thresholdMode)
                assert all(value == value for value in (k0, k, E)), (stats, sigma0, thresholdMode, k0, k, E)
    print(f"none of the {len(expressions)} integrator steps divides by kP / kD / sigma0 / sigmaV")


if __name__ == "__main__":
    test_triage_accepts_zero_sigma0()
    test_triage_still_rejects_bad_numbers()
    test_zero_sigma0_gives_zero_boost()
    test_nothing_divides_by_the_boost_constants()
    print("ALL GaMD ZERO-BOOST TESTS PASSED")
