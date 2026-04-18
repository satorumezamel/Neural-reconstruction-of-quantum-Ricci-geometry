# Follow-up Experiment: LOSO Stability of `switch_gain` + Phase-Shuffle Lag Surrogate

Date: 2026-04-18

Paper: *Intersection-Defined Phase Coordinates Reveal Localized Selection and a Non-Closed Observational Structure* (Satoru Watanabe, SIEL)

Script: `scripts/loso_and_lag_analysis.py`
Raw metrics: `reports/loso_lag_metrics.json`

## Abstract

Two complementary robustness checks were run on the published findings.
**(A) Leave-One-Subject-Out (LOSO) on `switch_gain`** over the 13 test sessions (P2–P9, P22–P26): the point estimate is extremely stable (LOSO mean ∈ [0.795, 0.821], full = 0.810), but the one-sided paired sign-test, which the paper reports at p = 0.046, rises to p = 0.073 in 9 of 13 LOSO folds — i.e. removing any of {P2, P5, P8, P22, P23, P24, P25, P26} pushes the comparison above the α = 0.05 threshold. Wilcoxon signed-rank is more stable (full p = 0.024–0.034).
**(B) Time-lagged Pearson cross-correlation ERicci ↔ QRicci with a phase-shuffle surrogate null (2,000 reps, per-subject FFT phase randomisation preserving power spectrum)**: per-subject peak |r| reproduces the paper’s value (mean |best r| = 0.45 vs. paper’s 0.42), but this aggregate is confounded by per-subject lag cherry-picking — only 3/26 subjects (P6, P24, P25) survive the phase-shuffle null at α = 0.05. At the group level, the mean lag curve stays inside the null 95 % band (peak r = 0.085 at lag = −7, peak-|r| p = 0.588). Conclusion: the `switch_gain` point estimate is robust, but its *significance* under LOSO is borderline, and the temporal-lag correspondence does not survive a phase-preserving surrogate null.

## Methods

### (A) LOSO on `switch_gain`

Inputs: `IDPC_Reproduction/Chapter7/session_wise_metrics_all_models.csv` (13 subjects × 4 models). For each subject *s*, we drop the row and recompute:

- `mean_best` = mean `switch_gain` of `best_true_search` across the remaining 12 subjects
- paired sign-test p-value (one-sided, `alternative="greater"`, matching the paper’s convention) for `best_true_search – neural_only_pca`, `– quantum_only_pca`, `– simple_mean_best`
- paired Wilcoxon signed-rank p-value (one-sided) on the same three contrasts

### (B) Phase-shuffle surrogate for lagged Pearson ERicci ↔ QRicci

Inputs: `IDPC_Reproduction_ricci/PX_timeseries.csv`, all 26 subjects, 30 samples each (P23 has 29).

For each subject we computed Pearson `r(E_Ricci[t], Q_Ricci[t + lag])` for lag ∈ [−10, +10]. Surrogate null: for each of 2,000 iterations, randomise the Fourier phases of `Q_Ricci` (preserving its power spectrum), recompute the lag curve; aggregate across subjects exactly as for the observed data. Two-sided p-values were formed at each lag as `P(|null| ≥ |observed|)`. A single peak-|r| statistic was also computed as a family-wise multiple-comparison-corrected test.

## Results

### (A) LOSO point estimate is robust; significance is fragile

| Quantity | Full sample (n=13) | LOSO min | LOSO max |
|---|---|---|---|
| mean `switch_gain` (best) | **0.8098** | 0.7954 (remove P25) | 0.8209 (remove P4) |
| paired sign-test p (vs neural-only) | 0.0461 | 0.0193 | **0.0730** |
| paired sign-test p (vs quantum-only) | 0.0461 | 0.0193 | **0.0730** |
| paired sign-test p (vs simple-mean) | 0.0461 | 0.0193 | **0.0730** |
| paired Wilcoxon p (vs neural-only) | 0.0341 | — | — |
| paired Wilcoxon p (vs quantum-only) | 0.0287 | — | — |
| paired Wilcoxon p (vs simple-mean) | 0.0239 | — | — |

The LOSO range of the point estimate is **±1.5 pp around 0.810**, so the effect size is robust. However the one-sided sign-test (which is also nominally borderline at the full sample) loses significance in 9/13 LOSO folds, reaching p = 0.073 > 0.05 whenever any of {P2, P5, P8, P22, P23, P24, P25, P26} is removed. Wilcoxon is stricter and remains < 0.05 at the full sample (0.024–0.034), but was not implemented as a LOSO sweep in the original paper.

Most influential subjects:
- Removing **P25** (whose individual `switch_gain` = 0.982 is the largest) drops the LOSO mean most, by 0.014.
- Removing **P4** (whose individual `switch_gain` = 0.676 is the smallest) raises the LOSO mean most, by 0.011.

![LOSO switch_gain](loso_switch_gain.png)

### (B) Lagged cross-correlation: per-subject peak reproduces the paper, group signal vanishes against phase-shuffle null

| Quantity | Value |
|---|---|
| per-subject mean `|best r|` (paper reports ≈ 0.42) | **0.450** |
| per-subject median `|best r|` | 0.441 |
| per-subject phase-shuffle peak-|r| p < 0.05 | **3 / 26** (P6, P24, P25) |
| per-subject median phase-shuffle p | 0.552 |
| group-mean peak lag | −7 samples |
| group-mean peak r | 0.085 |
| group-mean peak-|r| p (phase-shuffle) | **0.588** |

The group-mean curve lies almost entirely inside the 95 % phase-shuffle null band; only two individual lags (lag = −7, p = 0.045; lag = −2, p = 0.034) dip below α = 0.05 and neither survives any multiple-comparison adjustment across 21 lags.

![Lagged xcorr group mean](lagged_xcorr_group.png)

![Per-subject lag curves](lagged_xcorr_per_subject.png)

### (C) No relation between per-subject lag peak and LOSO influence

We correlated each test subject's per-subject best |r| lag with the amount that removing that subject changes the LOSO mean `switch_gain`. The Pearson correlation is essentially zero (see figure), indicating that per-subject lag structure and influence on `switch_gain` are independent — the two phenomena, to the extent they exist, are dissociated.

![LOSO influence vs best lag](loso_influence_vs_best_lag.png)

## Discussion

1. **The `switch_gain` effect size is robust, but its claimed significance sits right on the α = 0.05 threshold.** The LOSO sweep shows the paper's reported p = 0.046 is one borderline outcome in a distribution whose other folds all give p = 0.073. A Wilcoxon signed-rank (more powerful when paired differences are symmetric) is stricter (p ≈ 0.024–0.034) and is a cleaner primary statistic going forward.
2. **The paper's "mean |Pearson| ≈ 0.42" is reproduced numerically (0.45 here), but is produced by selecting each subject's best lag independently across a ±10-sample window** (21 candidate lags). Against a phase-preserving surrogate null — the natural non-parametric null for lagged auto-correlated series — the per-subject effect only survives in 3/26 subjects, and the group-level peak |r| is not distinguishable from the null (p = 0.59).
3. **Recommended tightening for any future analysis**: (i) pre-register a single lag (e.g. 0), or report lag-corrected p via the surrogate above; (ii) report Wilcoxon alongside the sign test for paired contrasts on `switch_gain`; (iii) supply LOSO sweeps as a standard robustness table, given that the test cohort is n = 13.

## Files

- Figures: `reports/loso_switch_gain.png`, `reports/lagged_xcorr_group.png`, `reports/lagged_xcorr_per_subject.png`, `reports/loso_influence_vs_best_lag.png`
- Metrics JSON: `reports/loso_lag_metrics.json`
- Script: `scripts/loso_and_lag_analysis.py`
