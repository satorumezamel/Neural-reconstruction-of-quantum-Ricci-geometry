"""
SIC (boundary-event) prediction experiment.

For each session:
  - Reconstruct the hybrid phi on the points axis (ported from
    IDPC_Reproduction.ipynb cell 16 / SECTION 1.8).
  - Compute features at each point t:
      phi(t), phi(t-1), phi(t-2), dphi(t),
      epsilon(t), depsilon(t),
      psi_E - psi_Q,  |psi_E - psi_Q|,
      time_since_last_event, last_event_J.
  - Build binary targets y_k(t) = 1 if any zero-crossing of h
    occurs in (t, t+k].  k in {1, 2, 3}.

Models (LOSO-CV across 26 sessions):
  - Logistic Regression
  - Random Forest
  - Baseline: time-since-last-event only (logistic).

Outputs:
  reports/sic_prediction_metrics.json
  reports/sic_prediction_auc.png
  reports/sic_prediction_features.png
  reports/sic_prediction_leadtime.png
"""
from __future__ import annotations

import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd

POINTS_DIR = "IDPC_Reproduction"
RICCI_DIR = "IDPC_Reproduction_ricci"
OUT_JSON = "reports/sic_prediction_metrics.json"
SESSIONS = [f"P{i}" for i in range(1, 27)]

# hybrid-phi params (kept identical to notebook cell 16)
ALPHA = 0.7
TAU = 3.0
KAPPA = 1.0
HYBRID_EPS = 1e-8
W_PRE = 1
W_POST = 2
EPS_SCALE = 0.5
EPS_MODE = "local"

LEAD_TIMES = [1, 2, 3]
N_RF = 200


# ---------- ported helpers from notebook cell 16 ----------
def finite_diff(x):
    x = np.asarray(x, float)
    dx = np.full_like(x, np.nan)
    n = len(x)
    if n < 2:
        return dx
    for i in range(1, n - 1):
        if np.isfinite(x[i - 1]) and np.isfinite(x[i + 1]):
            dx[i] = 0.5 * (x[i + 1] - x[i - 1])
    if n >= 2 and np.isfinite(x[0]) and np.isfinite(x[1]):
        dx[0] = x[1] - x[0]
    if n >= 2 and np.isfinite(x[-1]) and np.isfinite(x[-2]):
        dx[-1] = x[-1] - x[-2]
    return dx


def zero_crossings(h):
    h = np.asarray(h, float)
    m = np.isfinite(h)
    if m.sum() < 2:
        return np.array([], dtype=int)
    hh = h.copy()
    if not m.all():
        idx = np.arange(len(hh))
        hh[~m] = np.interp(idx[~m], idx[m], hh[m])
    s = np.sign(hh)
    return np.where((s[1:] * s[:-1]) < 0)[0]


def gaussian_delta_eps(h, eps):
    h = np.asarray(h, float)
    eps = float(max(eps, 1e-12))
    return np.exp(-(h ** 2) / (2 * eps ** 2)) / (np.sqrt(2 * np.pi) * eps)


def build_hybrid_phi(h, alpha=ALPHA, tau=TAU, kappa=KAPPA):
    h = np.asarray(h, float)
    n = len(h)
    phi_full = np.full(n, np.nan)
    J_tilde = np.full(n, np.nan)
    g_full = np.full(n, np.nan)
    J_event = np.full(n, np.nan)

    if n < 5:
        return phi_full, J_tilde, g_full, J_event, np.array([], dtype=int)

    dh = finite_diff(h)
    flips = zero_crossings(h)

    h_std_global = np.nanstd(h)
    eps_global = EPS_SCALE * h_std_global if np.isfinite(h_std_global) and h_std_global > 0 else 1e-6

    for k in flips:
        lo = max(0, k - W_PRE); hi = min(n, k + W_POST + 1)
        hw = h[lo:hi]; dhw = dh[lo:hi]
        if np.sum(np.isfinite(hw)) < 3:
            continue
        loc_std = np.nanstd(hw)
        eps_loc = EPS_SCALE * loc_std if (np.isfinite(loc_std) and loc_std > 0) else eps_global
        delta_eps = gaussian_delta_eps(hw, eps_loc)
        J_hat = np.nansum(dhw * delta_eps)
        J_event[k] = J_hat

    t_idx = np.arange(n, dtype=float)
    event_idx = np.where(np.isfinite(J_event))[0]
    if len(event_idx) == 0:
        J_tilde[:] = 0.0
    else:
        num = np.zeros(n); den = np.zeros(n)
        for k in event_idx:
            w = np.exp(-((t_idx - k) ** 2) / (2 * tau ** 2))
            num += J_event[k] * w
            den += w
        J_tilde = num / (den + HYBRID_EPS)

    s_h = np.nanstd(h); s_h = 1.0 if (not np.isfinite(s_h) or s_h < 1e-12) else s_h
    s_J = np.nanstd(J_tilde[np.isfinite(J_tilde)])
    s_J = 1.0 if (not np.isfinite(s_J) or s_J < 1e-12) else s_J
    rho = np.nanstd(np.abs(h))
    rho = 1.0 if (not np.isfinite(rho) or rho < 1e-12) else rho
    g_full = np.exp(-np.abs(h) / rho)

    base0 = (1.0 - g_full[0]) * (h[0] / s_h) + g_full[0] * kappa * (J_tilde[0] / s_J)
    phi_full[0] = base0
    for i in range(1, n):
        base = (1.0 - g_full[i]) * (h[i] / s_h) + g_full[i] * kappa * (J_tilde[i] / s_J)
        phi_full[i] = alpha * phi_full[i - 1] + (1.0 - alpha) * base
    return phi_full, J_tilde, g_full, J_event, flips


# ---------- per-session feature construction ----------
def wrap_pi(a):
    return np.angle(np.exp(1j * a))


def build_session_features(label):
    pts = pd.read_csv(os.path.join(POINTS_DIR, f"{label}_co_recon_features_W30_points.csv"))
    pts = pts[pts["valid"] == 1].reset_index(drop=True).copy()
    h = pts["h"].to_numpy(float)
    phi, J_tilde, g, J_event, flips = build_hybrid_phi(h)

    # ---- Ricci: map task_idx-axis E/Q onto points axis ----
    ricci_path = os.path.join(RICCI_DIR, f"{label}_timeseries.csv")
    eeg_path = os.path.join(POINTS_DIR, f"{label}_eeg_timeseries.csv")
    eeg = pd.read_csv(eeg_path, usecols=["task_idx"])
    ricci = pd.read_csv(ricci_path)

    # take the first len(h) rows of eeg task_idx (points == eeg rows that fall
    # into the valid window)
    t_eeg = eeg["task_idx"].astype(int).to_numpy()
    n_pts_in_eeg = min(len(t_eeg), len(h))
    E_map = dict(zip(ricci["task_idx"].astype(int), ricci["E_Ricci"].astype(float)))
    Q_map = dict(zip(ricci["task_idx"].astype(int), ricci["Q_Ricci"].astype(float)))
    E_pts = np.array([E_map.get(int(t), np.nan) for t in t_eeg[:n_pts_in_eeg]])
    Q_pts = np.array([Q_map.get(int(t), np.nan) for t in t_eeg[:n_pts_in_eeg]])

    # pad/truncate to match phi length
    if len(E_pts) < len(h):
        E_pts = np.concatenate([E_pts, np.full(len(h) - len(E_pts), np.nan)])
        Q_pts = np.concatenate([Q_pts, np.full(len(h) - len(Q_pts), np.nan)])
    else:
        E_pts = E_pts[: len(h)]
        Q_pts = Q_pts[: len(h)]

    # phase via Hilbert-like construction (x, dx/dt) plane
    dE = finite_diff(E_pts); dQ = finite_diff(Q_pts)
    psiE = np.arctan2(dE, E_pts); psiQ = np.arctan2(dQ, Q_pts)
    phase_diff = wrap_pi(psiE - psiQ)
    abs_phase_diff = np.abs(phase_diff)

    # epsilon(t) = phase_diff - smooth_phase_diff (structural residual)
    kern = np.ones(11) / 11.0
    smooth_pd = np.convolve(np.nan_to_num(phase_diff), kern, mode="same")
    eps = phase_diff - smooth_pd
    deps = np.concatenate([[np.nan], np.diff(eps)])

    dphi = np.concatenate([[np.nan], np.diff(phi)])
    phi_m1 = np.concatenate([[np.nan], phi[:-1]])
    phi_m2 = np.concatenate([[np.nan, np.nan], phi[:-2]])

    # time since last flip + last J
    n = len(h)
    last_event_dist = np.full(n, np.nan)
    last_J = np.full(n, 0.0)
    last_t = -1; last_val = 0.0
    flip_set = set(flips.tolist())
    for t in range(n):
        if last_t >= 0:
            last_event_dist[t] = t - last_t
        last_J[t] = last_val
        if t in flip_set:
            last_t = t
            last_val = float(J_event[t]) if np.isfinite(J_event[t]) else 0.0

    features = pd.DataFrame({
        "label": label,
        "t": np.arange(n),
        "h": h,
        "phi": phi,
        "phi_m1": phi_m1,
        "phi_m2": phi_m2,
        "dphi": dphi,
        "eps": eps,
        "deps": deps,
        "phase_diff": phase_diff,
        "abs_phase_diff": abs_phase_diff,
        "time_since_event": last_event_dist,
        "last_event_J": last_J,
    })
    # build targets at lead times 1..3 (event in (t, t+k])
    for k in LEAD_TIMES:
        y = np.zeros(n, dtype=int)
        for flip in flips:
            lo = flip - k; hi = flip
            lo = max(0, lo); hi = min(n, hi)
            if lo <= hi:
                y[lo:hi] = 1
        features[f"y_k{k}"] = y

    return features, flips


def build_all_features():
    out = []
    totals = {"n_points": 0, "n_events": 0}
    for lab in SESSIONS:
        try:
            df, flips = build_session_features(lab)
        except FileNotFoundError:
            continue
        totals["n_points"] += len(df)
        totals["n_events"] += len(flips)
        out.append(df)
    full = pd.concat(out, ignore_index=True)
    return full, totals


# ---------- modelling ----------
FEATURES_FULL = [
    "phi", "phi_m1", "phi_m2", "dphi", "eps", "deps",
    "phase_diff", "abs_phase_diff", "time_since_event", "last_event_J",
]
FEATURES_BASELINE = ["time_since_event"]


def fit_eval(X_tr, y_tr, X_te, y_te, model):
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score

    if model == "lr":
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    elif model == "rf":
        clf = RandomForestClassifier(n_estimators=N_RF, max_depth=12,
                                     class_weight="balanced_subsample",
                                     n_jobs=-1, random_state=42)
    elif model == "baseline":
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    else:
        raise ValueError(model)

    m_tr = np.isfinite(X_tr).all(axis=1); m_te = np.isfinite(X_te).all(axis=1)
    if m_tr.sum() < 50 or m_te.sum() < 50 or y_tr[m_tr].sum() == 0 or y_te[m_te].sum() == 0:
        return None
    clf.fit(X_tr[m_tr], y_tr[m_tr])
    proba = clf.predict_proba(X_te[m_te])[:, 1]
    pred = (proba >= 0.5).astype(int)
    y_te_m = y_te[m_te]
    return {
        "auc": float(roc_auc_score(y_te_m, proba)),
        "precision": float(precision_score(y_te_m, pred, zero_division=0)),
        "recall": float(recall_score(y_te_m, pred, zero_division=0)),
        "f1": float(f1_score(y_te_m, pred, zero_division=0)),
        "n_test": int(m_te.sum()),
        "n_pos_test": int(y_te_m.sum()),
        "feature_importance": (clf.feature_importances_.tolist()
                               if hasattr(clf, "feature_importances_") else None),
        "proba": proba.tolist(),
        "y_true": y_te_m.tolist(),
    }


def run_loso(full_df, target_col, feature_cols, model_name):
    sessions = sorted(full_df["label"].unique(),
                      key=lambda s: int(s[1:]))
    rows = []; importances = []
    for held_out in sessions:
        tr = full_df[full_df["label"] != held_out]
        te = full_df[full_df["label"] == held_out]
        X_tr = tr[feature_cols].to_numpy(float); y_tr = tr[target_col].to_numpy(int)
        X_te = te[feature_cols].to_numpy(float); y_te = te[target_col].to_numpy(int)
        r = fit_eval(X_tr, y_tr, X_te, y_te, model_name)
        if r is None:
            continue
        rows.append({
            "label": held_out, "auc": r["auc"], "precision": r["precision"],
            "recall": r["recall"], "f1": r["f1"],
            "n_test": r["n_test"], "n_pos_test": r["n_pos_test"],
        })
        if r["feature_importance"] is not None:
            importances.append(r["feature_importance"])
    rows_df = pd.DataFrame(rows)
    agg = {
        "mean_auc": float(rows_df["auc"].mean()) if len(rows_df) else float("nan"),
        "median_auc": float(rows_df["auc"].median()) if len(rows_df) else float("nan"),
        "std_auc": float(rows_df["auc"].std()) if len(rows_df) else float("nan"),
        "mean_precision": float(rows_df["precision"].mean()) if len(rows_df) else float("nan"),
        "mean_recall": float(rows_df["recall"].mean()) if len(rows_df) else float("nan"),
        "mean_f1": float(rows_df["f1"].mean()) if len(rows_df) else float("nan"),
        "n_sessions_scored": int(len(rows_df)),
    }
    if importances:
        mean_imp = np.mean(importances, axis=0)
        agg["feature_importance_mean"] = {
            f: float(v) for f, v in zip(feature_cols, mean_imp)
        }
    return agg, rows_df


# ---------- main ----------
def main():
    os.makedirs("reports", exist_ok=True)
    print("building features …")
    full, totals = build_all_features()
    print(f"  n_points={totals['n_points']}  n_events={totals['n_events']}")

    metrics = {"n_points_total": totals["n_points"],
               "n_events_total": totals["n_events"],
               "per_leadtime": {}}
    for k in LEAD_TIMES:
        metrics["per_leadtime"][f"k{k}"] = {}
        for model_name, feats in (("lr", FEATURES_FULL),
                                  ("rf", FEATURES_FULL),
                                  ("baseline", FEATURES_BASELINE)):
            agg, rows_df = run_loso(full, f"y_k{k}", feats, model_name)
            metrics["per_leadtime"][f"k{k}"][model_name] = agg
            rows_df.to_csv(f"reports/sic_prediction_k{k}_{model_name}_per_session.csv",
                           index=False)
            print(f"  k={k} model={model_name}  AUC={agg['mean_auc']:.3f}  "
                  f"(precision={agg['mean_precision']:.3f}, recall={agg['mean_recall']:.3f})")

    # --- |dphi| analysis: split test points into "high |dphi|" quartile and
    # compare AUC on each subset ---
    print("|dphi| quartile analysis …")
    k = 2
    target = f"y_k{k}"
    full2 = full.dropna(subset=["dphi"]).copy()
    full2["abs_dphi"] = full2["dphi"].abs()
    q = full2["abs_dphi"].quantile([0.25, 0.5, 0.75]).to_list()
    quart_results = {"thresholds": {"q25": q[0], "q50": q[1], "q75": q[2]}}
    for model in ("lr", "rf"):
        agg_high, _ = run_loso(full2[full2["abs_dphi"] >= q[2]].copy(),
                               target, FEATURES_FULL, model)
        agg_low, _ = run_loso(full2[full2["abs_dphi"] <= q[0]].copy(),
                              target, FEATURES_FULL, model)
        quart_results[model] = {"high_quartile_auc": agg_high["mean_auc"],
                                "low_quartile_auc": agg_low["mean_auc"]}
        print(f"  {model}: high-|dphi| AUC={agg_high['mean_auc']:.3f}, "
              f"low-|dphi| AUC={agg_low['mean_auc']:.3f}")
    metrics["abs_dphi_analysis"] = quart_results

    with open(OUT_JSON, "w") as f:
        json.dump(metrics, f, indent=2, default=float)
    print("wrote", OUT_JSON)

    # expose the feature matrix for downstream figures
    full.to_parquet("reports/sic_prediction_features.parquet",
                    compression="snappy") if False else None
    full.to_csv("reports/sic_prediction_features.csv", index=False)


if __name__ == "__main__":
    main()
