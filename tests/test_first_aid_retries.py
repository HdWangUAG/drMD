"""
Checks for the firstAid retry decorator in src/Surgery/drFirstAid.py.
Run with:  python tests/test_first_aid_retries.py   (or pytest tests/)

firstAidMaxRetries counts the firstAid retries allowed *after* the first attempt, so the wrapped
simulation function must always be run at least once, and at most firstAidMaxRetries + 1 times.
With firstAidMaxRetries: 0 the simulation is run once and is not rescued - it used to not be run
at all, and the decorator then reported a hard coded "Particle coordinate NaN" for a simulation
that had never started.

No OpenMM simulation is run here: the wrapped function and the firstAid protocol are both stubs.
"""
import sys
import re
from os import path as p

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from openmm import OpenMMException

from Surgery import drFirstAid
from UtilitiesCloset import drSplash, drSplicer

REAL_ERROR = "Energy is NaN after 37 steps of 09_cmd_stats"
FIRST_AID_SOURCE = p.join(ROOT, "src", "Surgery", "drFirstAid.py")


def make_kwargs(maxRetries):
    """The keyword arguments drOperator passes to a simulation function"""
    return {"prmtop": None,
            "inpcrd": None,
            "sim": {"stepName": "09_cmd_stats", "nSteps": 100},
            "saveFile": "start.xml",
            "outDir": p.join(ROOT, "tests", "no_such_out_dir"),
            "platform": None,
            "refPdb": "ref.pdb",
            "config": {"miscInfo": {"firstAidMaxRetries": maxRetries}}}


def stub_simulation(nFailures, error=None):
    """A simulation function that crashes nFailures times, then succeeds. Counts its own calls"""
    calls = []
    def simulationFunction(**kwargs):
        calls.append(kwargs["sim"]["stepName"])
        if len(calls) <= nFailures:
            raise error if error is not None else OpenMMException(REAL_ERROR)
        return "final.xml"
    return simulationFunction, calls


class Patched:
    """
    Replaces the parts of the firstAid protocol that would touch OpenMM or the file system:
    the firstAid protocol itself, the output merger and the splash screen.
    """
    def __enter__(self):
        self.firstAidCalls = []
        self.mergeCalls = []
        self.splashMessages = []
        self.originals = (drFirstAid.run_first_aid_protocol,
                          drSplicer.merge_partial_outputs,
                          drSplash.print_first_aid_failed,
                          drSplash.print_performing_first_aid)
        def fake_protocol(retries, maxRetries, *args, **kwargs):
            self.firstAidCalls.append(retries)
            return kwargs["saveFile"], retries + 1
        drFirstAid.run_first_aid_protocol = fake_protocol
        drSplicer.merge_partial_outputs = lambda **kwargs: self.mergeCalls.append(kwargs["simDir"])
        drSplash.print_first_aid_failed = lambda error: self.splashMessages.append(str(error))
        drSplash.print_performing_first_aid = lambda: None
        return self

    def __exit__(self, *exceptionInfo):
        (drFirstAid.run_first_aid_protocol,
         drSplicer.merge_partial_outputs,
         drSplash.print_first_aid_failed,
         drSplash.print_performing_first_aid) = self.originals
        return False


def run_wrapped(simulationFunction, maxRetries):
    """Decorate simulationFunction and call it, returning (result, raisedError)"""
    wrapped = drFirstAid.firstAid_handler()(simulationFunction)
    try:
        return wrapped(**make_kwargs(maxRetries)), None
    except Exception as error:
        return None, error


def test_single_handler_definition():
    """firstAid_handler was defined twice, the second definition silently won"""
    with open(FIRST_AID_SOURCE, "r") as fileHandle:
        definitions = re.findall(r"^def firstAid_handler", fileHandle.read(), flags=re.MULTILINE)
    print(f"firstAid_handler is defined {len(definitions)} time(s)")
    assert len(definitions) == 1, definitions


def test_zero_retries_runs_the_simulation_once():
    simulationFunction, calls = stub_simulation(nFailures=0)
    with Patched() as patched:
        result, error = run_wrapped(simulationFunction, maxRetries=0)
    print(f"firstAidMaxRetries=0, simulation succeeds: {len(calls)} call(s), returned {result}")
    assert error is None, error
    assert len(calls) == 1, calls
    assert result == "final.xml", result
    assert patched.firstAidCalls == [], "no firstAid should be attempted when the simulation works"


def test_zero_retries_does_not_retry_and_reports_the_real_error():
    simulationFunction, calls = stub_simulation(nFailures=99)
    with Patched() as patched:
        result, error = run_wrapped(simulationFunction, maxRetries=0)
    print(f"firstAidMaxRetries=0, simulation crashes: {len(calls)} call(s), raised {type(error).__name__}: {error}")
    assert len(calls) == 1, f"the simulation must run exactly once, ran {len(calls)} times"
    assert patched.firstAidCalls == [], "firstAidMaxRetries=0 means no firstAid retries"
    ## raising a string used to turn the real error into "exceptions must derive from BaseException"
    assert not isinstance(error, TypeError), error
    assert isinstance(error, OpenMMException), type(error)
    assert REAL_ERROR in str(error), error
    ## and the user must be shown the real reason, not a hard coded one
    assert patched.splashMessages == [REAL_ERROR], patched.splashMessages
    assert "Particle coordinate NaN" not in "".join(patched.splashMessages)


def test_zero_retries_reports_a_value_error_too():
    simulationFunction, calls = stub_simulation(nFailures=99, error=ValueError(REAL_ERROR))
    with Patched() as patched:
        result, error = run_wrapped(simulationFunction, maxRetries=0)
    print(f"firstAidMaxRetries=0, ValueError: raised {type(error).__name__}: {error}")
    assert len(calls) == 1, calls
    assert isinstance(error, ValueError) and REAL_ERROR in str(error), error


def test_retries_until_the_simulation_recovers():
    simulationFunction, calls = stub_simulation(nFailures=2)
    with Patched() as patched:
        result, error = run_wrapped(simulationFunction, maxRetries=3)
    print(f"firstAidMaxRetries=3, succeeds on attempt 3: {len(calls)} call(s), {len(patched.firstAidCalls)} firstAid attempt(s)")
    assert error is None, error
    assert result == "final.xml", result
    assert len(calls) == 3, calls
    assert patched.firstAidCalls == [0, 1], patched.firstAidCalls
    ## partial trajectories of the crashed attempts get merged on success
    assert len(patched.mergeCalls) == 1, patched.mergeCalls


def test_n_retries_run_the_simulation_n_plus_one_times():
    for maxRetries in [1, 2, 5]:
        simulationFunction, calls = stub_simulation(nFailures=99)
        with Patched() as patched:
            result, error = run_wrapped(simulationFunction, maxRetries=maxRetries)
        print(f"firstAidMaxRetries={maxRetries}, never recovers: {len(calls)} call(s), raised {type(error).__name__}")
        assert len(calls) == maxRetries + 1, f"{maxRetries} retries should give {maxRetries + 1} calls, got {len(calls)}"
        assert len(patched.firstAidCalls) == maxRetries, patched.firstAidCalls
        assert isinstance(error, OpenMMException) and REAL_ERROR in str(error), error


def test_unexpected_errors_are_not_retried():
    simulationFunction, calls = stub_simulation(nFailures=99, error=KeyError("phaseOfTheMoon"))
    with Patched() as patched:
        result, error = run_wrapped(simulationFunction, maxRetries=5)
    print(f"unexpected {type(error).__name__} raised straight through after {len(calls)} call(s)")
    assert isinstance(error, KeyError), error
    assert len(calls) == 1, calls
    assert patched.firstAidCalls == [], patched.firstAidCalls


if __name__ == "__main__":
    test_single_handler_definition()
    test_zero_retries_runs_the_simulation_once()
    test_zero_retries_does_not_retry_and_reports_the_real_error()
    test_zero_retries_reports_a_value_error_too()
    test_retries_until_the_simulation_recovers()
    test_n_retries_run_the_simulation_n_plus_one_times()
    test_unexpected_errors_are_not_retried()
    print("ALL FIRST AID TESTS PASSED")
