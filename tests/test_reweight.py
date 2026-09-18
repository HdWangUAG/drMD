"""
Checks for src/ExaminationRoom/drReweight.py on synthetic data with a known answer.
Run with:  python tests/test_reweight.py   (or pytest tests/)
"""
import sys, math
from os import path as p
import numpy as np

ROOT = p.dirname(p.dirname(p.abspath(__file__)))
sys.path.insert(0, p.join(ROOT, "src"))
from ExaminationRoom import drReweight

T = 300.0
kT = drReweight.KB_KCAL * T
beta = 1 / kT


def synthetic_gamd_sample(nFrames, seed):
    """
    Frames from a boosted ensemble along x with a known free energy F(x) = 2 (x^2 - 1)^2 kcal/mol.
    In the boosted ensemble the boost given x is Gaussian, dV | x ~ N(m(x), s^2). Like a real GaMD
    boost, m(x) is largest at the bottom of the wells, m(x) = 0.7 (8 - F(x)), so the boosted marginal
        p_b(x) ~ exp(-beta F(x) - beta m(x))  =  exp(-0.3 beta F(x)) * const
    is a flattened landscape, and second-order cumulant reweighting must recover F(x) up to a constant.
    """
    rng = np.random.default_rng(seed)
    F = lambda x: 2.0 * (x ** 2 - 1) ** 2
    m = lambda x: 0.7 * (8.0 - F(x))
    s = 1.5
    grid = np.linspace(-2, 2, 20001)
    logp = -beta * F(grid) - beta * m(grid)
    pdf = np.exp(logp - logp.max()); cdf = np.cumsum(pdf); cdf /= cdf[-1]
    x = np.interp(rng.random(nFrames), cdf, grid)
    dV = rng.normal(m(x), s)
    return x, dV, F


def test_cumulant_recovers_known_pmf():
    x, dV, F = synthetic_gamd_sample(400000, 1)
    result = drReweight.reweight_pmf(x[:, None], dV, T, bins=[40], minValues=[-1.6], maxValues=[1.6], method="cumulant2")
    centers = result["centers"][0]
    expected = F(centers); expected -= expected.min()
    defined = ~np.isnan(result["pmf"])
    err = np.abs(result["pmf"][defined] - expected[defined])
    print(f"cumulant2: max |F_est - F| = {err.max():.3f} kcal/mol over {defined.sum()} bins, anharmonicity = {result['anharmonicity']:.4f}")
    assert err.max() < 0.3, err.max()
    assert result["anharmonicity"] < drReweight.ANHARMONICITY_THRESHOLD
    ## without reweighting the profile is wrong by the boost term
    raw = drReweight.reweight_pmf(x[:, None], dV, T, bins=[40], minValues=[-1.6], maxValues=[1.6], method="none")
    rawErr = np.abs(raw["pmf"][defined] - expected[defined]).max()
    print(f"no reweighting: max |F_raw - F| = {rawErr:.3f} kcal/mol (should be large)")
    assert rawErr > 2.0
    ## exponential averaging is unbiased in the limit but noisier; just check it is in the right ballpark
    expResult = drReweight.reweight_pmf(x[:, None], dV, T, bins=[40], minValues=[-1.6], maxValues=[1.6], method="exp")
    expErr = np.abs(expResult["pmf"][defined] - expected[defined]).max()
    print(f"exp: max |F_est - F| = {expErr:.3f} kcal/mol (noisier than cumulant2)")
    assert expErr < 3.0


def test_minima_and_replicates():
    pmfs = []
    for seed in (2, 3):
        x, dV, F = synthetic_gamd_sample(200000, seed)
        pmfs.append(drReweight.reweight_pmf(x[:, None], dV, T, bins=[40], minValues=[-1.6], maxValues=[1.6]))
    minima = drReweight.find_minima(pmfs[0]["pmf"], pmfs[0]["centers"], maxMinima=2)
    positions = sorted(mn["position"][0] for mn in minima)
    print(f"minima found at {positions} (expected near -1 and +1)")
    assert len(positions) == 2 and abs(positions[0] + 1) < 0.15 and abs(positions[1] - 1) < 0.15
    comparison = drReweight.compare_replicates([pmf["pmf"] for pmf in pmfs])
    print(f"replicates: RMSD {comparison['rmsd']:.3f}, max |diff| {comparison['maxAbsDiff']:.3f} kcal/mol")
    assert comparison["rmsd"] < 0.2


def test_2d_periodic_grid():
    rng = np.random.default_rng(4)
    n = 400000
    phi = rng.uniform(-180, 180, n); psi = rng.uniform(-180, 180, n)
    dV = rng.normal(5, 1, n)
    result = drReweight.reweight_pmf(np.column_stack([phi, psi]), dV, T, bins=[12, 12], minValues=[-180, -180], maxValues=[180, 180])
    assert result["pmf"].shape == (12, 12)
    assert result["counts"].sum() == n, "every frame must land in exactly one bin, including the +180 edge"
    assert np.nanmax(result["pmf"]) < 0.35, "uniform sampling with x-independent boost must give a flat PMF"
    print("2D grid OK: flat PMF within", round(float(np.nanmax(result["pmf"])), 3), "kcal/mol")


if __name__ == "__main__":
    test_cumulant_recovers_known_pmf()
    test_minima_and_replicates()
    test_2d_periodic_grid()
    print("ALL REWEIGHT TESTS PASSED")
