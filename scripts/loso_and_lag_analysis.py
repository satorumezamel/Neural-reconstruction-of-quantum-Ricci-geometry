"""
Follow-up experiment combining two robustness analyses:

  (A) Leave-One-Subject-Out (LOSO) stability of `switch_gain`
      -- using session_wise_metrics_all_models.csv from Chapter 7 --

  (B) Time-lagged Pearson cross-correlation ERicci -> QRicci
      with a phase-shuffle surrogate null
      -- using IDPC_Reproduction_ricci/PX_timeseries.csv --

Outputs (under reports/):
  - loso_lag_metrics.json
  - loso_switch_gain.png
  - loso_influence_vs_best_lag.png
  - lagged_xcorr_group.png
  - lagged_xcorr_per_subject.png
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
METRICS_CSV = ROOT / "IDPC_Reproduction" / "Chapter7" / "session_wise_metrics_all_models.csv"
RICCI_DIR = ROOT / "IDPC_Reproduction_ricci"
REPORTS = ROOT / "reports"
REPORTS.mkdir(exist_ok=True)

RNG = np.random.default_rng(2026_04_18)


# ---------- (A) LOSO on switch_gain -----------------------------------------

def sign_test_p(x: np.ndarray, alternative: str = "greater") -> float:
    """Sign test p-value.  Default `greater` matches the paper's convention
    (one-sided: best-model differences are positive)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x != 0)]
    if len(x) == 0:
        return np.nan
    k = int((x > 0).sum())
    n = len(x)
    return float(stats.binomtest(k, n, 0.5, alternative=alternative).pvalue)


def wilcoxon_p(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x != 0)]
    if len(x) < 2:
        return np.nan
    try:
        return float(stats.wilcoxon(x, alternative="greater").pvalue)
    except ValueError:
        return np.nan


def loso_switch_gain():
    df = pd.read_csv(METRICS_CSV)
    df["subject"] = df["label"].str.extract(r"(P\d+)", expand=False)

    models = ["best_true_search", "neural_only_pca", "quantum_only_pca", "simple_mean_best"]
    piv = df.pivot_table(index="subject", columns="model", values="switch_gain")
    piv = piv[models].dropna()
    subjects = piv.index.tolist()

    # full sample baseline
    full_mean = {m: float(piv[m].mean()) for m in models}
    full_paired = {
        "best_minus_neural":  sign_test_p(piv["best_true_search"] - piv["neural_only_pca"]),
        "best_minus_quantum": sign_test_p(piv["best_true_search"] - piv["quantum_only_pca"]),
        "best_minus_simple":  sign_test_p(piv["best_true_search"] - piv["simple_mean_best"]),
        "wilcoxon_best_minus_neural":  wilcoxon_p(piv["best_true_search"] - piv["neural_only_pca"]),
        "wilcoxon_best_minus_quantum": wilcoxon_p(piv["best_true_search"] - piv["quantum_only_pca"]),
        "wilcoxon_best_minus_simple":  wilcoxon_p(piv["best_true_search"] - piv["simple_mean_best"]),
    }

    # LOSO
    records = []
    for s in subjects:
        sub = piv.drop(index=s)
        rec = {
            "removed": s,
            "n_remaining": int(len(sub)),
            "mean_best": float(sub["best_true_search"].mean()),
            "mean_neural": float(sub["neural_only_pca"].mean()),
            "mean_quantum": float(sub["quantum_only_pca"].mean()),
            "mean_simple": float(sub["simple_mean_best"].mean()),
            "p_best_vs_neural": sign_test_p(sub["best_true_search"] - sub["neural_only_pca"]),
            "p_best_vs_quantum": sign_test_p(sub["best_true_search"] - sub["quantum_only_pca"]),
            "p_best_vs_simple":  sign_test_p(sub["best_true_search"] - sub["simple_mean_best"]),
        }
        rec["influence_on_best"] = float(full_mean["best_true_search"] - rec["mean_best"])
        records.append(rec)
    loso_df = pd.DataFrame(records).sort_values("mean_best").reset_index(drop=True)

    return {
        "subjects": subjects,
        "piv": piv,
        "full_mean": full_mean,
        "full_paired": full_paired,
        "loso_df": loso_df,
    }


# ---------- (B) Lagged Pearson + phase-shuffle surrogate --------------------

def phase_shuffle(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """FFT phase-randomisation preserving the power spectrum."""
    x = np.asarray(x, float)
    n = len(x)
    X = np.fft.rfft(x - x.mean())
    amp = np.abs(X)
    # randomise phases of all but DC; at Nyquist keep real when n even
    phases = rng.uniform(0.0, 2 * np.pi, size=len(X))
    phases[0] = 0.0
    if n % 2 == 0:
        phases[-1] = 0.0
    Xs = amp * np.exp(1j * phases)
    y = np.fft.irfft(Xs, n=n)
    return y + x.mean()


def pearson_at_lag(a: np.ndarray, b: np.ndarray, lag: int) -> float:
    """Correlation between a[t] and b[t+lag] (positive lag = b leads)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if lag >= 0:
        x = a[: len(a) - lag]; y = b[lag:]
    else:
        x = a[-lag:]; y = b[: len(b) + lag]
    if len(x) < 5 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def lagged_xcorr(max_lag: int = 10, n_surrogate: int = 2000):
    files = sorted(RICCI_DIR.glob("P*_timeseries.csv"),
                   key=lambda p: int(re.search(r"P(\d+)", p.stem).group(1)))
    lags = np.arange(-max_lag, max_lag + 1)

    per_subj = {}
    per_subj_curve = []
    null_mean_curves = []

    for fp in files:
        subj = re.search(r"(P\d+)", fp.stem).group(1)
        d = pd.read_csv(fp).dropna(subset=["E_Ricci", "Q_Ricci"])
        E = d["E_Ricci"].to_numpy(float)
        Q = d["Q_Ricci"].to_numpy(float)
        if len(E) < max_lag * 2 + 5:
            continue
        curve = np.array([pearson_at_lag(E, Q, int(l)) for l in lags])
        per_subj_curve.append(curve)

        # per-subject surrogate: phase-shuffle Q
        null = np.zeros((n_surrogate, len(lags)))
        for k in range(n_surrogate):
            Qs = phase_shuffle(Q, RNG)
            null[k] = [pearson_at_lag(E, Qs, int(l)) for l in lags]
        per_subj[subj] = {
            "curve": curve.tolist(),
            "best_lag_abs": int(lags[np.nanargmax(np.abs(curve))]),
            "best_r_abs":   float(curve[np.nanargmax(np.abs(curve))]),
            "null_max_abs_95": float(np.nanpercentile(np.nanmax(np.abs(null), axis=1), 95)),
            "p_peak_abs": float((np.nanmax(np.abs(null), axis=1) >= np.nanmax(np.abs(curve))).mean()),
        }
        null_mean_curves.append(null.mean(axis=0))

    per_subj_curve = np.vstack(per_subj_curve)
    group = per_subj_curve.mean(axis=0)

    # group-level null: average null curves the same way (per-subject shuffle is
    # independent, so mean-of-means is the correct aggregate)
    # For a stronger group null we re-shuffle all subjects jointly:
    n_subj = per_subj_curve.shape[0]
    group_null = np.zeros((n_surrogate, len(lags)))
    # Re-read once to save I/O
    cached = []
    for fp in files:
        d = pd.read_csv(fp).dropna(subset=["E_Ricci", "Q_Ricci"])
        cached.append((d["E_Ricci"].to_numpy(float), d["Q_Ricci"].to_numpy(float)))
    cached = cached[: n_subj]
    for k in range(n_surrogate):
        curves = []
        for (E, Q) in cached:
            Qs = phase_shuffle(Q, RNG)
            curves.append([pearson_at_lag(E, Qs, int(l)) for l in lags])
        group_null[k] = np.nanmean(curves, axis=0)

    # two-sided p at each lag
    p_two_sided = ((np.abs(group_null) >= np.abs(group)).sum(axis=0) + 1) / (n_surrogate + 1)
    p_peak_abs = ((np.nanmax(np.abs(group_null), axis=1) >= np.nanmax(np.abs(group))) .sum() + 1) / (n_surrogate + 1)

    return {
        "lags": lags.tolist(),
        "group_curve": group.tolist(),
        "group_null_mean": group_null.mean(axis=0).tolist(),
        "group_null_lo": np.percentile(group_null, 2.5, axis=0).tolist(),
        "group_null_hi": np.percentile(group_null, 97.5, axis=0).tolist(),
        "p_two_sided": p_two_sided.tolist(),
        "p_peak_abs": float(p_peak_abs),
        "per_subj": per_subj,
        "per_subj_curves": per_subj_curve.tolist(),
        "subjects": [re.search(r"(P\d+)", fp.stem).group(1) for fp in files][: n_subj],
    }


# ---------- Plots ------------------------------------------------------------

def plot_loso(loso):
    piv = loso["piv"]
    loso_df = loso["loso_df"]
    full_mean = loso["full_mean"]

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    # Left: LOSO mean switch_gain (best_true_search) per removed subject
    ax[0].axhline(full_mean["best_true_search"], color="k", ls="--",
                  label=f"full-sample mean = {full_mean['best_true_search']:.3f}")
    ax[0].bar(loso_df["removed"], loso_df["mean_best"], color="#3a7ca5")
    ax[0].set_ylabel("LOSO mean switch_gain (best_true_search)")
    ax[0].set_xlabel("removed subject")
    ax[0].set_title("LOSO stability of mean switch_gain")
    ax[0].tick_params(axis="x", rotation=45)
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    # Right: paired sign-test p-values vs each baseline across LOSO
    for col, lbl, c in [
        ("p_best_vs_neural",  "best vs neural-only",  "#e07a5f"),
        ("p_best_vs_quantum", "best vs quantum-only", "#81b29a"),
        ("p_best_vs_simple",  "best vs simple-mean",  "#f2cc8f"),
    ]:
        ax[1].plot(loso_df["removed"], loso_df[col], "o-", label=lbl, color=c)
    ax[1].axhline(0.05, color="k", ls=":", label="alpha = 0.05")
    ax[1].set_ylabel("paired sign-test p (LOSO)")
    ax[1].set_xlabel("removed subject")
    ax[1].set_title("Stability of paired comparisons")
    ax[1].tick_params(axis="x", rotation=45)
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    fig.tight_layout()
    out = REPORTS / "loso_switch_gain.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return str(out)


def plot_influence_vs_lag(loso, lagres):
    loso_df = loso["loso_df"].set_index("removed")
    per = lagres["per_subj"]

    xs, ys, labels = [], [], []
    for subj, row in loso_df.iterrows():
        if subj in per:
            xs.append(per[subj]["best_lag_abs"])
            ys.append(row["influence_on_best"])
            labels.append(subj)
    xs = np.array(xs); ys = np.array(ys)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.axhline(0, color="k", lw=0.8)
    ax.axvline(0, color="k", lw=0.8)
    ax.scatter(xs, ys, color="#3d405b", s=45)
    for x, y, lb in zip(xs, ys, labels):
        ax.annotate(lb, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=8)
    if len(xs) > 3:
        r, p = stats.pearsonr(xs, ys)
        ax.set_title(f"LOSO influence on switch_gain vs per-subject best lag\n"
                     f"Pearson r = {r:.2f}, p = {p:.3f}")
    ax.set_xlabel("per-subject best |r| lag (samples)")
    ax.set_ylabel("LOSO influence on switch_gain\n(full - LOSO)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = REPORTS / "loso_influence_vs_best_lag.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return str(out)


def plot_lag_curves(lagres):
    lags = np.array(lagres["lags"])
    group = np.array(lagres["group_curve"])
    lo = np.array(lagres["group_null_lo"])
    hi = np.array(lagres["group_null_hi"])
    p = np.array(lagres["p_two_sided"])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(lags, lo, hi, alpha=0.25, color="#888", label="phase-shuffle null 95% band")
    ax.axhline(0, color="k", lw=0.8)
    ax.plot(lags, group, "o-", color="#3d405b", lw=2, label="observed group-mean Pearson r")
    for lg, rv, pv in zip(lags, group, p):
        if pv < 0.05:
            ax.plot(lg, rv, "*", color="#e63946", ms=12)
    ax.set_xlabel("lag (samples, positive = Q leads)")
    ax.set_ylabel("group-mean Pearson r(E_Ricci(t), Q_Ricci(t+lag))")
    ax.set_title(f"Time-lagged correlation across {len(lagres['subjects'])} subjects\n"
                 f"phase-shuffle peak-|r| p = {lagres['p_peak_abs']:.3f}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = REPORTS / "lagged_xcorr_group.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return str(out)


def plot_per_subject_lag(lagres):
    curves = np.array(lagres["per_subj_curves"])
    lags = np.array(lagres["lags"])
    subjects = lagres["subjects"]

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, s in enumerate(subjects):
        ax.plot(lags, curves[i], lw=0.8, alpha=0.6)
    ax.plot(lags, curves.mean(0), "k-", lw=2.5, label="mean")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("lag (samples)")
    ax.set_ylabel("Pearson r(E_Ricci(t), Q_Ricci(t+lag))")
    ax.set_title("Per-subject lag curves (thin) + group mean (thick)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = REPORTS / "lagged_xcorr_per_subject.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return str(out)


# ---------- Main -------------------------------------------------------------

def main():
    print("[1/2] LOSO switch_gain analysis ...")
    loso = loso_switch_gain()
    print(f"      subjects: {len(loso['subjects'])}, full mean best = {loso['full_mean']['best_true_search']:.4f}")

    print("[2/2] Lagged cross-correlation with phase-shuffle surrogate ...")
    lagres = lagged_xcorr(max_lag=10, n_surrogate=2000)
    lags = np.array(lagres["lags"])
    peak_idx = int(np.nanargmax(np.abs(lagres["group_curve"])))
    print(f"      peak lag = {lags[peak_idx]}, r = {lagres['group_curve'][peak_idx]:.3f}, "
          f"peak-|r| p = {lagres['p_peak_abs']:.4f}")

    f1 = plot_loso(loso)
    f2 = plot_influence_vs_lag(loso, lagres)
    f3 = plot_lag_curves(lagres)
    f4 = plot_per_subject_lag(lagres)
    print("plots:", f1, f2, f3, f4, sep="\n  ")

    # serialisable metrics
    metrics = {
        "loso": {
            "subjects": loso["subjects"],
            "full_mean": loso["full_mean"],
            "full_paired_p": loso["full_paired"],
            "loso_records": loso["loso_df"].to_dict(orient="records"),
            "range_mean_best": [float(loso["loso_df"]["mean_best"].min()),
                                float(loso["loso_df"]["mean_best"].max())],
            "max_p_best_vs_neural": float(loso["loso_df"]["p_best_vs_neural"].max()),
            "max_p_best_vs_quantum": float(loso["loso_df"]["p_best_vs_quantum"].max()),
            "max_p_best_vs_simple":  float(loso["loso_df"]["p_best_vs_simple"].max()),
            "most_influential_positive": loso["loso_df"].iloc[0]["removed"],   # smallest LOSO mean (largest positive influence)
            "most_influential_negative": loso["loso_df"].iloc[-1]["removed"],  # largest LOSO mean (negative influence)
        },
        "lag": {
            "lags": lagres["lags"],
            "group_curve": lagres["group_curve"],
            "group_null_lo": lagres["group_null_lo"],
            "group_null_hi": lagres["group_null_hi"],
            "p_two_sided_per_lag": lagres["p_two_sided"],
            "peak_lag": int(lags[peak_idx]),
            "peak_r": float(lagres["group_curve"][peak_idx]),
            "peak_abs_p": float(lagres["p_peak_abs"]),
            "n_subjects": len(lagres["subjects"]),
            "per_subj_best_lag": {k: v["best_lag_abs"] for k, v in lagres["per_subj"].items()},
            "per_subj_best_r":   {k: v["best_r_abs"]   for k, v in lagres["per_subj"].items()},
            "per_subj_p":        {k: v["p_peak_abs"]   for k, v in lagres["per_subj"].items()},
        },
    }
    out_json = REPORTS / "loso_lag_metrics.json"
    with open(out_json, "w") as f:
        json.dump(metrics, f, indent=2)
    print("metrics saved to", out_json)


if __name__ == "__main__":
    main()
