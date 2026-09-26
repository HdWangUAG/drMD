## BASIC PYTHON LIBRARIES
import os
from os import path as p
import argparse
import json
import math
import warnings

## NUMERICAL LIBRARIES
import numpy as np
import pandas as pd

## CLEAN CODE
from typing import List, Optional, Tuple, Dict, Sequence

########################################################################################################
## Boltzmann constant in kcal/(mol K)
KB_KCAL: float = 0.0019872041
## anharmonicity below which the second-order cumulant expansion is considered reliable (Miao et al. 2014)
ANHARMONICITY_THRESHOLD: float = 0.01
## frames per bin below which the cumulant term -0.5 * beta * var(dV) is too noisy to trust, and the
## default minCount. Measured on alanine dipeptide, two independent GaMD replicates, 36x36 bins,
## cumulant2: replicate-to-replicate PMF RMSD / max |diff| in kcal/mol was 2.26 / 6.77 at 10 frames
## per bin, 0.58 / 1.70 at 100, 0.41 / 1.17 at 300 and 0.22 / 0.74 at 1000.
MIN_FRAMES_PER_BIN: int = 300
########################################################################################################
"""
Reweighting of Gaussian accelerated MD (GaMD) simulations run with drMD.

A GaMD trajectory samples a boosted potential V' = V + dV. The unbiased distribution along any
coordinate A is recovered by reweighting each frame with exp(beta * dV):

    p*(A) = p(A) * <exp(beta dV)>_A / <exp(beta dV)>

Because dV is a sum of many contributions its distribution is close to Gaussian, and the
exponential average is far better estimated by a cumulant expansion truncated at second order
(Miao, Sinko, Pierce, Bucher, Walker & McCammon, JCTC 2014, 10, 2677):

    ln <exp(beta dV)>_A  ~  beta <dV>_A + (beta^2 / 2) var(dV)_A

so that the free energy along A is

    F(A) = -kT ln p(A) - <dV>_A - (beta / 2) var(dV)_A + const

Direct exponential averaging is implemented for reference only: at the boost magnitudes used in
practice it is dominated by the few largest dV values and is very noisy.

Three diagnostics decide whether a reweighted profile is usable and are reported with every PMF:
  * the anharmonicity of the dV distribution, gamma = S_max - S (0 for a Gaussian); the cumulant
    expansion is reliable for gamma below ~0.01
  * the occupancy of the bins actually used: how many bins were kept, how many dropped, and the
    smallest and median number of frames among the kept bins. The variance of dV per bin is the
    noisiest ingredient of the expansion, so bins holding fewer than MIN_FRAMES_PER_BIN frames are
    dropped by default and a thinner profile is warned about
  * agreement between independent replicates

Inputs are the gamd.log written by drGaMD (columns step, time_ps, V_total, V_dihedral, dV_P, dV_D in
kcal/mol) and a cv.csv written by drCVReporter (columns step, time_ps, cv_0, ...), which have one row
per saved frame each.

Command line usage (one PMF per replicate, plus their comparison):

    python drReweight.py --gamdLog rep1/gamd.log rep2/gamd.log \\
                         --cv rep1/00_reporters_and_plots/cv.csv rep2/00_reporters_and_plots/cv.csv \\
                         --columns cv_0 cv_1 --temperature 300 --bins 36 --min -180 -180 --max 180 180 \\
                         --periodic --outDir reweighted
"""
########################################################################################################
def read_gamd_log(gamdLog: str) -> pd.DataFrame:
    """Reads a gamd.log written by drGaMD. Adds a dV column = dV_P + dV_D (kcal/mol)."""
    gamdDf: pd.DataFrame = pd.read_csv(gamdLog)
    gamdDf["dV"] = gamdDf["dV_P"] + gamdDf["dV_D"]
    return gamdDf
########################################################################################################
def read_cv_csv(cvCsv: str) -> pd.DataFrame:
    """Reads a cv.csv written by drCVReporter."""
    return pd.read_csv(cvCsv)
########################################################################################################
def align_frames(gamdDf: pd.DataFrame, cvDf: pd.DataFrame) -> pd.DataFrame:
    """
    Joins gamd.log and cv.csv on the step column. Both are written at logInterval, so normally the
    join is one-to-one; frames present in only one file (e.g. after a resume) are dropped.
    Duplicate steps (a resumed run re-reporting a frame) keep the last occurrence.
    """
    gamdDf = gamdDf.drop_duplicates(subset="step", keep="last")
    cvDf = cvDf.drop_duplicates(subset="step", keep="last")
    merged: pd.DataFrame = pd.merge(gamdDf, cvDf.drop(columns=[c for c in ["time_ps"] if c in cvDf.columns]),
                                    on="step", how="inner")
    return merged
########################################################################################################
def anharmonicity(dV: np.ndarray, nBins: int = 50) -> float:
    """
    Anharmonicity of the boost potential distribution, gamma = S_max - S (Miao et al. 2014), where
    S is the differential entropy of the dV distribution estimated from a histogram and S_max the
    entropy of a Gaussian with the same variance. gamma is 0 for a Gaussian and grows as the
    distribution becomes skewed or heavy-tailed.

    Args:
        dV (np.ndarray): boost potentials (kcal/mol), one per frame
        nBins (int): number of histogram bins

    Returns:
        gamma (float)
    """
    dV = np.asarray(dV, dtype=float)
    sigma: float = dV.std()
    if sigma == 0.0 or len(dV) < 2:
        return 0.0
    counts, edges = np.histogram(dV, bins=nBins, density=True)
    binWidth: float = edges[1] - edges[0]
    nonZero = counts[counts > 0]
    entropy: float = -np.sum(nonZero * np.log(nonZero) * binWidth)
    entropyMax: float = 0.5 * math.log(2 * math.pi * math.e * sigma ** 2)
    return float(entropyMax - entropy)
########################################################################################################
def make_bin_edges(cv: np.ndarray,
                    bins: Sequence[int],
                      minValues: Optional[Sequence[float]] = None,
                        maxValues: Optional[Sequence[float]] = None) -> List[np.ndarray]:
    """Bin edges along each coordinate, defaulting to the sampled range."""
    nDims: int = cv.shape[1]
    edges: List[np.ndarray] = []
    for dim in range(nDims):
        lo: float = cv[:, dim].min() if minValues is None else minValues[dim]
        hi: float = cv[:, dim].max() if maxValues is None else maxValues[dim]
        edges.append(np.linspace(lo, hi, bins[dim] + 1))
    return edges
########################################################################################################
def reweight_pmf(cv: np.ndarray,
                  dV: np.ndarray,
                    temperature: float,
                      bins: Sequence[int],
                        minValues: Optional[Sequence[float]] = None,
                          maxValues: Optional[Sequence[float]] = None,
                            method: str = "cumulant2",
                              minCount: int = MIN_FRAMES_PER_BIN) -> Dict:
    """
    Potential of mean force along one or more coordinates from a GaMD trajectory.

    Args:
        cv (np.ndarray): shape (nFrames, nDims) coordinate values per frame
        dV (np.ndarray): shape (nFrames,) total boost potential per frame (kcal/mol)
        temperature (float): simulation temperature (K)
        bins (Sequence[int]): number of bins along each coordinate
        minValues, maxValues: histogram range per coordinate (default: sampled range)
        method (str): "cumulant2" (default, second-order cumulant expansion), "exp" (direct exponential
                      average, reference only) or "none" (no reweighting: the boosted distribution)
        minCount (int): bins with fewer frames than this are left undefined (NaN). The default,
                        MIN_FRAMES_PER_BIN, is the point at which replicates agree; lowering it warns.

    Returns:
        dict with keys:
            pmf (np.ndarray): free energy in kcal/mol on the bin grid, shifted so its minimum is 0
            edges (list of np.ndarray): bin edges per coordinate
            centers (list of np.ndarray): bin centres per coordinate
            counts (np.ndarray): frames per bin
            meanDV (np.ndarray), varDV (np.ndarray): per-bin boost statistics (kcal/mol, (kcal/mol)^2)
            anharmonicity (float): gamma of the whole dV distribution
            minCount (int): the threshold used
            nBinsKept (int), nBinsDropped (int): bins defined and left undefined
            minCountKept (int), medianCountKept (float): frames in the thinnest and the median kept bin
            occupancyOk (bool): every kept bin holds at least MIN_FRAMES_PER_BIN frames
            method (str), temperature (float), nFrames (int)
    """
    cv = np.atleast_2d(np.asarray(cv, dtype=float))
    if cv.shape[0] == 1 and cv.shape[1] > 1 and len(dV) == cv.shape[1]:
        cv = cv.T
    dV = np.asarray(dV, dtype=float)
    if len(dV) != cv.shape[0]:
        raise ValueError(f"cv has {cv.shape[0]} frames but dV has {len(dV)}")
    nDims: int = cv.shape[1]
    if len(bins) != nDims:
        raise ValueError(f"bins must have one entry per coordinate ({nDims}), got {len(bins)}")
    kT: float = KB_KCAL * temperature
    beta: float = 1.0 / kT

    edges: List[np.ndarray] = make_bin_edges(cv, bins, minValues, maxValues)
    ## bin index of every frame along every coordinate; frames outside the range are dropped
    binIndexes: np.ndarray = np.zeros((cv.shape[0], nDims), dtype=int)
    inRange: np.ndarray = np.ones(cv.shape[0], dtype=bool)
    for dim in range(nDims):
        idx: np.ndarray = np.digitize(cv[:, dim], edges[dim]) - 1
        ## the top edge belongs to the last bin
        idx[cv[:, dim] == edges[dim][-1]] = bins[dim] - 1
        inRange &= (idx >= 0) & (idx < bins[dim])
        binIndexes[:, dim] = idx
    binIndexes = binIndexes[inRange]
    dVUsed: np.ndarray = dV[inRange]
    flatIndex: np.ndarray = np.ravel_multi_index(binIndexes.T, tuple(bins))
    nBinsTotal: int = int(np.prod(bins))

    counts: np.ndarray = np.bincount(flatIndex, minlength=nBinsTotal).astype(float)
    sumDV: np.ndarray = np.bincount(flatIndex, weights=dVUsed, minlength=nBinsTotal)
    sumDV2: np.ndarray = np.bincount(flatIndex, weights=dVUsed ** 2, minlength=nBinsTotal)
    with np.errstate(divide="ignore", invalid="ignore"):
        meanDV: np.ndarray = sumDV / counts
        varDV: np.ndarray = sumDV2 / counts - meanDV ** 2
        probability: np.ndarray = counts / counts.sum()
        if method == "cumulant2":
            ## F = -kT [ln p + beta <dV> + beta^2 var/2]
            pmf: np.ndarray = -kT * np.log(probability) - meanDV - 0.5 * beta * varDV
        elif method == "exp":
            ## F = -kT ln( sum_i exp(beta dV_i) / N ), computed with a log-sum-exp for stability
            shift: float = beta * dVUsed.max() if len(dVUsed) > 0 else 0.0
            expSum: np.ndarray = np.bincount(flatIndex, weights=np.exp(beta * dVUsed - shift), minlength=nBinsTotal)
            pmf: np.ndarray = -kT * (np.log(expSum) + shift - math.log(len(dVUsed)))
        elif method == "none":
            pmf: np.ndarray = -kT * np.log(probability)
        else:
            raise ValueError(f"method must be 'cumulant2', 'exp' or 'none', got {method}")
    kept: np.ndarray = counts >= minCount
    pmf[~kept] = np.nan
    if np.all(np.isnan(pmf)):
        raise ValueError(f"No bin holds the {minCount} frames needed to define a free energy; "
                         "use fewer bins, run longer, or lower minCount and accept a noisier profile")
    pmf -= np.nanmin(pmf)

    ## per-PMF occupancy, so a profile built from nearly empty bins is visible rather than inferred
    keptCounts: np.ndarray = counts[kept]
    minCountKept: int = int(keptCounts.min())
    occupancyOk: bool = bool(minCountKept >= MIN_FRAMES_PER_BIN)
    if not occupancyOk:
        warnings.warn(f"thin PMF: the sparsest of the {int(kept.sum())} bins used holds only {minCountKept} frames "
                      f"(minCount = {minCount}, {MIN_FRAMES_PER_BIN} recommended). The cumulant term needs the "
                      "variance of dV per bin, so this profile is dominated by noise even if it looks smooth.",
                      stacklevel=2)

    shape: Tuple[int, ...] = tuple(bins)
    centers: List[np.ndarray] = [0.5 * (edge[1:] + edge[:-1]) for edge in edges]
    return {"pmf": pmf.reshape(shape),
            "edges": edges,
            "centers": centers,
            "counts": counts.reshape(shape),
            "meanDV": meanDV.reshape(shape),
            "varDV": varDV.reshape(shape),
            "anharmonicity": anharmonicity(dVUsed) if len(dVUsed) > 1 else 0.0,
            "minCount": int(minCount),
            "nBinsKept": int(kept.sum()),
            "nBinsDropped": int(nBinsTotal - kept.sum()),
            "minCountKept": minCountKept,
            "medianCountKept": float(np.median(keptCounts)),
            "occupancyOk": occupancyOk,
            "method": method,
            "temperature": temperature,
            "nFrames": int(len(dVUsed))}
########################################################################################################
def compare_replicates(pmfs: List[np.ndarray]) -> Dict:
    """
    Agreement between PMFs from independent replicates computed on the same grid.
    Each PMF is shifted to a minimum of 0 first. Bins undefined (NaN) in any replicate are ignored.

    Returns:
        dict with mean (np.ndarray), std (np.ndarray, standard deviation across replicates per bin),
        rmsd (float, root mean square deviation of every replicate from the mean profile, kcal/mol),
        maxAbsDiff (float, largest deviation of any replicate from the mean in any bin, kcal/mol),
        nReplicates (int), nBinsCompared (int)
    """
    if len(pmfs) < 2:
        raise ValueError("compare_replicates needs at least two PMFs")
    stack: np.ndarray = np.stack([pmf - np.nanmin(pmf) for pmf in pmfs])
    defined: np.ndarray = ~np.any(np.isnan(stack), axis=0)
    mean: np.ndarray = np.where(defined, np.nanmean(stack, axis=0), np.nan)
    std: np.ndarray = np.where(defined, np.nanstd(stack, axis=0), np.nan)
    deviations: np.ndarray = (stack - mean)[:, defined]
    return {"mean": mean,
            "std": std,
            "rmsd": float(np.sqrt(np.mean(deviations ** 2))) if deviations.size else float("nan"),
            "maxAbsDiff": float(np.abs(deviations).max()) if deviations.size else float("nan"),
            "nReplicates": len(pmfs),
            "nBinsCompared": int(defined.sum())}
########################################################################################################
def find_minima(pmf: np.ndarray, centers: List[np.ndarray], maxMinima: int = 5, minSeparation: int = 2) -> List[Dict]:
    """
    Local minima of a PMF (1D or 2D grid), lowest first. Two minima closer than minSeparation bins
    along every axis are reported once. Useful for checking minima positions against another method.
    """
    minima: List[Dict] = []
    finitePmf: np.ndarray = np.where(np.isnan(pmf), np.inf, pmf)
    order: np.ndarray = np.argsort(finitePmf, axis=None)
    for flat in order:
        if not np.isfinite(finitePmf.flat[flat]):
            break
        index: Tuple[int, ...] = np.unravel_index(flat, pmf.shape)
        ## must be no higher than every defined neighbour
        isMin: bool = True
        for dim in range(pmf.ndim):
            for delta in (-1, 1):
                neighbour = list(index)
                neighbour[dim] = (neighbour[dim] + delta) % pmf.shape[dim]
                if finitePmf[tuple(neighbour)] < finitePmf[index]:
                    isMin = False
        if not isMin:
            continue
        if any(all(abs(index[dim] - other["index"][dim]) < minSeparation for dim in range(pmf.ndim)) for other in minima):
            continue
        minima.append({"index": tuple(int(i) for i in index),
                       "position": [float(centers[dim][index[dim]]) for dim in range(pmf.ndim)],
                       "freeEnergy": float(pmf[index])})
        if len(minima) >= maxMinima:
            break
    return minima
########################################################################################################
def write_pmf_csv(result: Dict, outCsv: str) -> None:
    """Writes a PMF as a long-format CSV: one row per bin with coordinates, free energy, counts and dV statistics."""
    grids = np.meshgrid(*result["centers"], indexing="ij")
    df = pd.DataFrame({f"cv_{dim}": grid.ravel() for dim, grid in enumerate(grids)})
    df["freeEnergy_kcal"] = result["pmf"].ravel()
    df["counts"] = result["counts"].ravel()
    df["meanDV_kcal"] = result["meanDV"].ravel()
    df["varDV_kcal2"] = result["varDV"].ravel()
    df.to_csv(outCsv, index=False)
########################################################################################################
def plot_pmf(result: Dict, outPng: str, title: str = "", replicateStd: Optional[np.ndarray] = None) -> None:
    """Plots a 1D or 2D PMF to PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pmf: np.ndarray = result["pmf"]
    centers: List[np.ndarray] = result["centers"]
    fig, ax = plt.subplots(figsize=(6, 5))
    if pmf.ndim == 1:
        ax.plot(centers[0], pmf, marker="o", ms=3, color="#1f77b4")
        if replicateStd is not None:
            ax.fill_between(centers[0], pmf - replicateStd, pmf + replicateStd, alpha=0.3, color="#1f77b4", label="replicate s.d.")
            ax.legend()
        ax.set_xlabel("cv_0"); ax.set_ylabel("free energy (kcal/mol)")
    elif pmf.ndim == 2:
        mesh = ax.pcolormesh(centers[0], centers[1], pmf.T, shading="nearest", cmap="viridis")
        fig.colorbar(mesh, ax=ax, label="free energy (kcal/mol)")
        ax.set_xlabel("cv_0"); ax.set_ylabel("cv_1")
    else:
        plt.close(fig)
        return
    ax.set_title(title or f"{result['method']} reweighted PMF (gamma = {result['anharmonicity']:.4f})")
    fig.tight_layout(); fig.savefig(outPng, dpi=150); plt.close(fig)
########################################################################################################
def reweight_protocol(gamdLogs: List[str],
                       cvCsvs: List[str],
                         columns: List[str],
                           temperature: float,
                             bins: Sequence[int],
                               minValues: Optional[Sequence[float]],
                                 maxValues: Optional[Sequence[float]],
                                   method: str,
                                     outDir: str,
                                       minCount: int = MIN_FRAMES_PER_BIN) -> Dict:
    """
    Reweights each replicate (one gamd.log + one cv.csv), writes PMF CSVs and plots, and compares
    the replicates. Returns a summary dictionary (also written to reweight_summary.json).
    """
    if len(gamdLogs) != len(cvCsvs):
        raise ValueError("Give one cv.csv per gamd.log")
    os.makedirs(outDir, exist_ok=True)
    results: List[Dict] = []
    summary: Dict = {"method": method, "temperature": temperature, "columns": columns, "bins": list(bins),
                     "anharmonicityThreshold": ANHARMONICITY_THRESHOLD,
                     "minCount": minCount, "minFramesPerBin": MIN_FRAMES_PER_BIN, "replicates": []}
    for repIndex, (gamdLog, cvCsv) in enumerate(zip(gamdLogs, cvCsvs)):
        merged: pd.DataFrame = align_frames(read_gamd_log(gamdLog), read_cv_csv(cvCsv))
        cv: np.ndarray = merged[columns].to_numpy()
        result: Dict = reweight_pmf(cv, merged["dV"].to_numpy(), temperature, bins, minValues, maxValues, method, minCount)
        results.append(result)
        tag: str = f"rep{repIndex + 1}"
        write_pmf_csv(result, p.join(outDir, f"pmf_{tag}.csv"))
        plot_pmf(result, p.join(outDir, f"pmf_{tag}.png"), title=f"{tag}: {method} (gamma = {result['anharmonicity']:.4f})")
        minima: List[Dict] = find_minima(result["pmf"], result["centers"])
        summary["replicates"].append({"gamdLog": gamdLog, "cvCsv": cvCsv, "nFrames": result["nFrames"],
                                      "anharmonicity": result["anharmonicity"],
                                      "anharmonicityOk": result["anharmonicity"] < ANHARMONICITY_THRESHOLD,
                                      "nBinsKept": result["nBinsKept"], "nBinsDropped": result["nBinsDropped"],
                                      "minCountKept": result["minCountKept"], "medianCountKept": result["medianCountKept"],
                                      "occupancyOk": result["occupancyOk"],
                                      "meanDV_kcal": float(merged["dV"].mean()), "stdDV_kcal": float(merged["dV"].std()),
                                      "minima": minima})
    if len(results) > 1:
        comparison: Dict = compare_replicates([result["pmf"] for result in results])
        summary["replicateAgreement"] = {"rmsd_kcal": comparison["rmsd"], "maxAbsDiff_kcal": comparison["maxAbsDiff"],
                                         "nBinsCompared": comparison["nBinsCompared"]}
        meanResult: Dict = dict(results[0]); meanResult["pmf"] = comparison["mean"]
        write_pmf_csv(meanResult, p.join(outDir, "pmf_mean.csv"))
        plot_pmf(meanResult, p.join(outDir, "pmf_mean.png"), title=f"mean of {len(results)} replicates ({method})",
                 replicateStd=comparison["std"] if comparison["mean"].ndim == 1 else None)
    with open(p.join(outDir, "reweight_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return summary
########################################################################################################
def main() -> None:
    parser = argparse.ArgumentParser(description="Reweight drMD GaMD simulations to a potential of mean force.")
    parser.add_argument("--gamdLog", nargs="+", required=True, help="gamd.log file(s), one per replicate")
    parser.add_argument("--cv", nargs="+", required=True, help="cv.csv file(s), one per replicate, same order")
    parser.add_argument("--columns", nargs="+", default=["cv_0"], help="cv.csv columns to build the PMF along")
    parser.add_argument("--temperature", type=float, default=300.0, help="temperature (K)")
    parser.add_argument("--bins", nargs="+", type=int, default=None, help="bins per coordinate (default 36 each)")
    parser.add_argument("--min", nargs="+", type=float, default=None, help="lower histogram edge per coordinate")
    parser.add_argument("--max", nargs="+", type=float, default=None, help="upper histogram edge per coordinate")
    parser.add_argument("--method", default="cumulant2", choices=["cumulant2", "exp", "none"])
    parser.add_argument("--minCount", type=int, default=MIN_FRAMES_PER_BIN,
                        help="bins with fewer frames are left undefined (default: %(default)s, below which "
                             "the cumulant term is too noisy for replicates to agree)")
    parser.add_argument("--outDir", default="reweighted")
    args = parser.parse_args()
    bins = args.bins if args.bins is not None else [36] * len(args.columns)
    summary = reweight_protocol(args.gamdLog, args.cv, args.columns, args.temperature, bins, args.min, args.max,
                                args.method, args.outDir, args.minCount)
    for rep in summary["replicates"]:
        flag = "OK" if rep["anharmonicityOk"] else "WARNING: dV distribution is not Gaussian enough for cumulant reweighting"
        print(f"{rep['gamdLog']}: {rep['nFrames']} frames, <dV> = {rep['meanDV_kcal']:.2f} +/- {rep['stdDV_kcal']:.2f} kcal/mol, "
              f"anharmonicity = {rep['anharmonicity']:.4f} [{flag}]")
        occupancyFlag = "OK" if rep["occupancyOk"] else f"WARNING: bins this thin do not reproduce between replicates, raise --minCount to {MIN_FRAMES_PER_BIN} or use fewer bins"
        print(f"    bins: {rep['nBinsKept']} kept, {rep['nBinsDropped']} dropped below minCount = {summary['minCount']}; "
              f"frames per kept bin: min {rep['minCountKept']}, median {rep['medianCountKept']:.0f} [{occupancyFlag}]")
        for minimum in rep["minima"]:
            print(f"    minimum at {minimum['position']} : {minimum['freeEnergy']:.2f} kcal/mol")
    if "replicateAgreement" in summary:
        agreement = summary["replicateAgreement"]
        print(f"replicate agreement: RMSD {agreement['rmsd_kcal']:.2f} kcal/mol, max |diff| {agreement['maxAbsDiff_kcal']:.2f} kcal/mol "
              f"over {agreement['nBinsCompared']} bins")
    print(f"outputs written to {args.outDir}")
########################################################################################################
if __name__ == "__main__":
    main()
