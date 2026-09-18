"""
Checks for the GAMD / metadynamics validation in src/Triage/drConfigTriage.py.
Run with:  python tests/test_config_triage.py   (or pytest tests/)
"""
import sys, copy, io, contextlib
from os import path as p

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from Triage import drConfigTriage

ATOM = {"CHAIN_ID": "A", "RES_NAME": "ALA", "RES_ID": 2, "ATOM_NAME": "CA"}
CUSTOM = {"keyword": "custom", "customSelection": [ATOM]}


def gamd_step(name, stage, **extra):
    info = {"stage": stage}
    info.update(extra)
    return {"stepName": name, "simulationType": "GAMD", "duration": "1 ns", "timestep": "2 fs", "temperature": 300,
            "logInterval": "10 ps", "gamdInfo": info}


def base_config(steps):
    return {"pathInfo": {"inputDir": ROOT, "outputDir": p.join(ROOT, "outputs")},
            "hardwareInfo": {"parallelCPU": 1, "platform": "CPU", "subprocessCpus": 1},
            "miscInfo": {"pH": 7, "firstAidMaxRetries": 1, "boxGeometry": "cubic", "writeMyMethodsSection": False,
                         "skipPdbTriage": True, "trajectorySelections": [{"selection": {"keyword": "all"}}]},
            "simulationInfo": steps}


def validate(config):
    """Returns the validated config, or None if validation failed (drSplash calls exit(1) after printing)."""
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            return drConfigTriage.validate_config(copy.deepcopy(config))
        except SystemExit:
            return None


def test_gamd_defaults_and_sequence():
    cfg = validate(base_config([gamd_step("01_stats", "cmd_stats"), gamd_step("02_equil", "gamd_equil"), gamd_step("03_prod", "gamd_prod")]))
    assert cfg is not None
    info = cfg["simulationInfo"][0]["gamdInfo"]
    assert info == {"stage": "cmd_stats", "boostType": "dual", "thresholdMode": "lower", "ensemble": "NPT",
                    "sigma0P": 6.0, "sigma0D": 6.0, "updateInterval": 500, "excludeRestraintsFromBoost": True, "cvs": []}, info
    print("GAMD defaults OK:", info)


def test_gamd_rejections():
    ## first GAMD step must be cmd_stats
    assert validate(base_config([gamd_step("01_equil", "gamd_equil")])) is None
    ## stage required / valid
    assert validate(base_config([gamd_step("01_x", "warmup")])) is None
    ## boost definition must not change between stages
    assert validate(base_config([gamd_step("01_stats", "cmd_stats", boostType="dual"), gamd_step("02_prod", "gamd_prod", boostType="total")])) is None
    ## ensemble must not change between stages
    assert validate(base_config([gamd_step("01_stats", "cmd_stats", ensemble="NVT"), gamd_step("02_prod", "gamd_prod")])) is None
    ## temperatureRange not allowed
    step = gamd_step("01_stats", "cmd_stats"); del step["temperature"]; step["temperatureRange"] = [300, 310]
    assert validate(base_config([step])) is None
    ## bad numbers
    assert validate(base_config([gamd_step("01_stats", "cmd_stats", sigma0P=-1)])) is None
    assert validate(base_config([gamd_step("01_stats", "cmd_stats", updateInterval=0)])) is None
    ## cvs use the bias syntax without grid keys; com_distance needs selection2
    ok = validate(base_config([gamd_step("01_stats", "cmd_stats", cvs=[{"biasVar": "com_distance", "selection": CUSTOM, "selection2": CUSTOM}])]))
    assert ok is not None
    assert validate(base_config([gamd_step("01_stats", "cmd_stats", cvs=[{"biasVar": "com_distance", "selection": CUSTOM}])])) is None
    ## unknown simulationType
    bad = gamd_step("01_stats", "cmd_stats"); bad["simulationType"] = "REMD"
    assert validate(base_config([bad])) is None
    print("GAMD rejections OK")


def meta_step(**extra):
    info = {"height": 0.8, "biasFactor": 10, "biases": [{"biasVar": "torsion", "minValue": -180, "maxValue": 180, "biasWidth": 5.73,
                                                          "selection": {"keyword": "custom", "customSelection": [ATOM] * 4}}]}
    info.update(extra)
    return {"stepName": "01_meta", "simulationType": "META", "duration": "1 ns", "timestep": "2 fs", "temperature": 300,
            "logInterval": "10 ps", "metaDynamicsInfo": info}


def test_metadynamics_frequency():
    cfg = validate(base_config([meta_step()]))
    assert cfg is not None
    info = cfg["simulationInfo"][0]["metaDynamicsInfo"]
    assert info["frequency"] == 500 and info["saveFrequency"] == 500, info
    cfg = validate(base_config([meta_step(frequency=250)]))
    assert cfg["simulationInfo"][0]["metaDynamicsInfo"]["saveFrequency"] == 250
    assert validate(base_config([meta_step(frequency=500, saveFrequency=1200)])) is None      # not a multiple
    assert validate(base_config([meta_step(frequency=500, saveFrequency=1000)])) is not None
    assert validate(base_config([meta_step(frequency=0)])) is None
    assert validate(base_config([meta_step(freeEnergyInterval="1 ns")])) is not None
    assert validate(base_config([meta_step(freeEnergyInterval="soon")])) is None
    ## the old bug: invalid metaDynamicsInfo (no biases) must fail validation
    bad = meta_step(); del bad["metaDynamicsInfo"]["biases"]
    assert validate(base_config([bad])) is None
    ## COM_DISTANCE bias with selection2
    com = meta_step(biases=[{"biasVar": "COM_DISTANCE", "minValue": 2, "maxValue": 10, "biasWidth": 0.5, "selection": CUSTOM, "selection2": CUSTOM}])
    assert validate(base_config([com])) is not None
    print("metadynamics frequency / saveFrequency / COM_DISTANCE validation OK")


if __name__ == "__main__":
    test_gamd_defaults_and_sequence()
    test_gamd_rejections()
    test_metadynamics_frequency()
    print("ALL CONFIG TRIAGE TESTS PASSED")
