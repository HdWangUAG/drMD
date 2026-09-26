"""
Checks that drMD.main reports a botched simulation in its exit status, not only on the terminal.

A batch where one run raises must exit non-zero and must NOT print "Simulations Complete!", so that
a wrapper script on a remote machine can tell a finished wave from a half-finished one. A batch where
every run succeeds must still exit 0 and still print it.

No simulation is run: drOperator.drMD_protocol is replaced, as are config reading and cleanup, and
main() is driven in a subprocess so its exit status can be read.

Run with:  python tests/test_exit_status.py   (or pytest tests/)
"""
import os
import sys
import shutil
import subprocess
import tempfile
from os import path as p

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
SRC = p.join(ROOT, "src")

SUCCESS_MESSAGE = "Simulations Complete!"

## driver run in a subprocess: mock out everything except drMD.main's own control flow
DRIVER = '''
import os, sys
sys.path.insert(0, {src!r})
import drMD
from Triage import drConfigTriage, drConfigWriter
from ExaminationRoom import drCleanup
from Surgery import drOperator

outDir, inDir, failingArg, nCpusArg = sys.argv[1:5]
failing = [name for name in failingArg.split(",") if name]
nCpus = int(nCpusArg)

batchConfig = {{"pathInfo": {{"inputDir": inDir, "outputDir": outDir}},
                "hardwareInfo": {{"parallelCPU": nCpus, "subprocessCpus": 1}},
                "miscInfo": {{"skipPdbTriage": True, "writeMyMethodsSection": False}}}}

drConfigTriage.read_input_yaml = lambda *args, **kwargs: batchConfig
drConfigTriage.validate_config = lambda config: config
drConfigWriter.make_per_protein_config = lambda pdbFile, config: pdbFile
drCleanup.clean_up_handler = lambda config: None
## keep the parallel branch in-process: only the wiring of the botched list is under test
if nCpus > 1:
    drMD.process_map = lambda function, items, max_workers=None: [function(item) for item in items]

def fake_protocol(runConfigYaml):
    name = os.path.splitext(os.path.basename(runConfigYaml))[0]
    if name in failing:
        raise TypeError("write_pdb() got an unexpected keyword argument (simulated crash)")

drOperator.drMD_protocol = fake_protocol

drMD.main("batchConfig.yaml")
'''


def run_batch(failing, nCpus=1):
    """Runs drMD.main over two mock systems, failing the named ones. Returns (returnCode, output)."""
    tmpDir = tempfile.mkdtemp(prefix="drmd_exit_status_")
    try:
        inDir = p.join(tmpDir, "inputs")
        outDir = p.join(tmpDir, "outputs")
        os.makedirs(inDir)
        os.makedirs(outDir)
        for name in ("systemA", "systemB"):
            with open(p.join(inDir, f"{name}.pdb"), "w") as pdbFile:
                pdbFile.write("END\n")
        driverPath = p.join(tmpDir, "driver.py")
        with open(driverPath, "w") as driver:
            driver.write(DRIVER.format(src=SRC))
        completed = subprocess.run([sys.executable, driverPath, outDir, inDir, ",".join(failing), str(nCpus)],
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=tmpDir)
        return completed.returncode, completed.stdout
    finally:
        shutil.rmtree(tmpDir, ignore_errors=True)


def test_botched_serial_run_exits_non_zero():
    returnCode, output = run_batch(["systemA"])
    print(f"serial, one botched run: exit code {returnCode}")
    assert returnCode != 0, f"a botched simulation must not exit 0\n{output}"
    assert SUCCESS_MESSAGE not in output, f"must not claim success when a run botched\n{output}"
    assert "systemA" in output, f"the failing system should be named\n{output}"


def test_successful_serial_run_exits_zero():
    returnCode, output = run_batch([])
    print(f"serial, no botched runs: exit code {returnCode}")
    assert returnCode == 0, f"a clean batch must exit 0\n{output}"
    assert SUCCESS_MESSAGE in output, f"a clean batch must still report success\n{output}"


def test_botched_parallel_run_exits_non_zero():
    returnCode, output = run_batch(["systemB"], nCpus=2)
    print(f"parallel, one botched run: exit code {returnCode}")
    assert returnCode != 0, f"a botched simulation must not exit 0 in parallel either\n{output}"
    assert SUCCESS_MESSAGE not in output, f"must not claim success when a run botched\n{output}"


if __name__ == "__main__":
    test_botched_serial_run_exits_non_zero()
    test_successful_serial_run_exits_zero()
    test_botched_parallel_run_exits_non_zero()
    print("ALL EXIT STATUS TESTS PASSED")
