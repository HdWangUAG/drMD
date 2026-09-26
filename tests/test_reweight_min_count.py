"""
Checks that drReweight does not quietly hand back a PMF built from nearly empty bins.

The cumulant term -0.5 * beta * var(dV) needs hundreds of frames per bin: measured on alanine
dipeptide with two independent GaMD replicates (36x36 bins, cumulant2), the replicate-to-replicate
PMF RMSD was 2.26 kcal/mol at minCount = 10 (max |diff| 6.77, one replicate putting its global
minimum in a sterically forbidden region) against 0.41 / 1.17 kcal/mol at minCount = 300.

So the default must not be 10, and a profile whose kept bins are still thinner than
drReweight.MIN_FRAMES_PER_BIN must be flagged, not returned quietly.

No MD is run: the dV and CV arrays are synthetic and the bin occupancy is exact.

Run with:  python tests/test_reweight_min_count.py   (or pytest tests/)
"""
import inspect
import re
import sys
import warnings
from os import path as p
import numpy as np

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from ExaminationRoom import drReweight

T = 300.0
REWEIGHT_SOURCE = p.join(ROOT, "src", "ExaminationRoom", "drReweight.py")


def occupancy_sample(perBinCounts, seed=0):
    """CV and dV arrays placing exactly perBinCounts[i] frames in the middle of bin i."""
    rng = np.random.default_rng(seed)
    cv = np.concatenate([np.full(count, index + 0.5) for index, count in enumerate(perBinCounts)])
    dV = rng.normal(8.0, 1.5, len(cv))
    return cv[:, None], dV


def reweight(perBinCounts, **kwargs):
    cv, dV = occupancy_sample(perBinCounts)
    return drReweight.reweight_pmf(cv, dV, T, bins=[len(perBinCounts)],
                                   minValues=[0], maxValues=[len(perBinCounts)], **kwargs)


def test_default_min_count_is_defensible():
    threshold = drReweight.MIN_FRAMES_PER_BIN
    print(f"drReweight.MIN_FRAMES_PER_BIN = {threshold}")
    assert threshold >= 300, f"the default must be high enough for the cumulant variance, got {threshold}"
    for function in (drReweight.reweight_pmf, drReweight.reweight_protocol):
        default = inspect.signature(function).parameters["minCount"].default
        print(f"{function.__name__} default minCount = {default}")
        assert default != 10, f"{function.__name__} still defaults to minCount = 10"
        assert default == threshold, f"{function.__name__} default should be MIN_FRAMES_PER_BIN, got {default}"
    ## the documented command line must not reintroduce the old default
    with open(REWEIGHT_SOURCE) as sourceFile:
        source = sourceFile.read()
    argparseLine = re.search(r'add_argument\("--minCount".*', source).group(0)
    print(f"argparse: {argparseLine.strip()}")
    assert "MIN_FRAMES_PER_BIN" in argparseLine, f"--minCount must default to MIN_FRAMES_PER_BIN: {argparseLine}"


def test_well_sampled_pmf_is_not_flagged():
    result = reweight([400] * 20)
    print(f"400 frames per bin: kept {result['nBinsKept']}, dropped {result['nBinsDropped']}, "
          f"min {result['minCountKept']}, median {result['medianCountKept']}, ok = {result['occupancyOk']}")
    assert result["occupancyOk"] is True
    assert result["nBinsKept"] == 20 and result["nBinsDropped"] == 0
    assert result["minCountKept"] == 400 and result["medianCountKept"] == 400


def test_dropped_bins_are_counted():
    ## half the bins are too thin for the default and must be reported as dropped, not just NaN
    result = reweight([400] * 10 + [50] * 10)
    print(f"half thin bins: kept {result['nBinsKept']}, dropped {result['nBinsDropped']}, "
          f"min {result['minCountKept']}, median {result['medianCountKept']}")
    assert result["nBinsKept"] == 10 and result["nBinsDropped"] == 10
    assert result["minCountKept"] == 400
    assert np.isnan(result["pmf"][15]), "a bin below minCount must stay undefined"


def test_thin_pmf_is_flagged():
    ## a user who lowers minCount below the threshold still gets a PMF, but must be told it is thin
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = reweight([25] * 20, minCount=20)
    print(f"25 frames per bin at minCount=20: ok = {result['occupancyOk']}, warnings = {[str(w.message) for w in caught]}")
    assert result["occupancyOk"] is False, "a PMF whose bins hold 25 frames must not be reported as ok"
    assert result["minCountKept"] == 25 and result["medianCountKept"] == 25
    assert result["nBinsKept"] == 20 and result["nBinsDropped"] == 0
    assert len(caught) > 0, "a thin PMF must raise a warning, not be returned quietly"
    assert any("minCount" in str(w.message) or "frames per bin" in str(w.message) for w in caught), \
        [str(w.message) for w in caught]


def test_thin_data_at_the_default_is_not_silent():
    ## with the default minCount there are no usable bins at all, which must be an error, not a garbage PMF
    try:
        reweight([25] * 20)
    except ValueError as error:
        print(f"thin data at the default raises: {error}")
    else:
        raise AssertionError("25 frames per bin must not silently give a PMF at the default minCount")


if __name__ == "__main__":
    test_default_min_count_is_defensible()
    test_well_sampled_pmf_is_not_flagged()
    test_dropped_bins_are_counted()
    test_thin_pmf_is_flagged()
    test_thin_data_at_the_default_is_not_silent()
    print("ALL REWEIGHT MIN COUNT TESTS PASSED")
