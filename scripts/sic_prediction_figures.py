"""Generate figures for the SIC-prediction experiment."""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPORTS = "reports"
METRICS = os.path.join(REPORTS, "sic_prediction_metrics.json")


def main():
    m = json.load(open(METRICS))

    # --- Figure 1: AUC by model, grouped over lead times ---
    lead_times = [1, 2, 3]
    models = ["baseline", "lr", "rf"]
    pretty = {"baseline": "Baseline\n(time-since-event)",
              "lr": "Logistic Regression",
              "rf": "Random Forest"}
    colors = {"baseline": "#b0b0b0", "lr": "#1f77b4", "rf": "#2ca02c"}
    x = np.arange(len(lead_times))
    w = 0.25
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, mdl in enumerate(models):
        means = [m["per_leadtime"][f"k{k}"][mdl]["mean_auc"] for k in lead_times]
        stds = [m["per_leadtime"][f"k{k}"][mdl]["std_auc"] for k in lead_times]
        ax.bar(x + (i - 1) * w, means, w, yerr=stds, capsize=3,
               color=colors[mdl], label=pretty[mdl])
    ax.set_xticks(x); ax.set_xticklabels([f"t+{k}" for k in lead_times])
    ax.set_ylabel("LOSO-CV mean AUC  (± SD across 26 sessions)")
    ax.set_ylim(0.45, 0.85)
    ax.axhline(0.5,  ls=":",  color="red", alpha=0.5, label="Chance (0.50)")
    ax.axhline(0.65, ls="--", color="red", alpha=0.5, label="H1 threshold (0.65)")
    ax.set_title("SIC prediction AUC by model and lead time")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    p1 = os.path.join(REPORTS, "sic_prediction_auc.png")
    plt.savefig(p1, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", p1)

    # --- Figure 2: Random-Forest feature importance at k=2 ---
    imp = m["per_leadtime"]["k2"]["rf"]["feature_importance_mean"]
    names = list(imp.keys()); vals = list(imp.values())
    order = np.argsort(vals)[::-1]
    names = [names[i] for i in order]
    vals = [vals[i] for i in order]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(range(len(names))[::-1], vals, color="#2ca02c", alpha=0.8)
    ax.set_yticks(range(len(names))[::-1])
    ax.set_yticklabels(names)
    ax.set_xlabel("Random Forest mean feature importance (k=t+2)")
    ax.set_title("Feature importance for SIC prediction  (RF, target = event within 2 points)")
    for i, v in enumerate(vals):
        ax.text(v + 0.003, len(names) - 1 - i, f"{v:.3f}",
                va="center", fontsize=9)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    p2 = os.path.join(REPORTS, "sic_prediction_features.png")
    plt.savefig(p2, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", p2)

    # --- Figure 3: lead time degradation & |dphi| panel ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # lead-time AUC curves per model
    for mdl in models:
        aucs = [m["per_leadtime"][f"k{k}"][mdl]["mean_auc"] for k in lead_times]
        stds = [m["per_leadtime"][f"k{k}"][mdl]["std_auc"] for k in lead_times]
        axes[0].errorbar(lead_times, aucs, yerr=stds, marker="o",
                         capsize=4, label=pretty[mdl], color=colors[mdl])
    axes[0].axhline(0.65, ls="--", color="red", alpha=0.5, label="H1 threshold (0.65)")
    axes[0].axhline(0.5,  ls=":",  color="red", alpha=0.5)
    axes[0].set_xticks(lead_times); axes[0].set_xlabel("Lead time (points ahead)")
    axes[0].set_ylabel("LOSO-CV mean AUC")
    axes[0].set_title("A. AUC vs lead time")
    axes[0].set_ylim(0.45, 0.85)
    axes[0].grid(True, alpha=0.3); axes[0].legend(loc="lower right", fontsize=9)

    # |dphi| quartile panel
    qa = m["abs_dphi_analysis"]
    xq = np.arange(2)
    lbls = ["low |dφ|\nquartile (≤q25)", "high |dφ|\nquartile (≥q75)"]
    w = 0.35
    for i, mdl in enumerate(("lr", "rf")):
        vals = [qa[mdl]["low_quartile_auc"], qa[mdl]["high_quartile_auc"]]
        axes[1].bar(xq + (i - 0.5) * w, vals, w, label=pretty[mdl],
                    color=colors[mdl])
    axes[1].axhline(0.65, ls="--", color="red", alpha=0.5, label="H1 threshold (0.65)")
    axes[1].axhline(0.5,  ls=":",  color="red", alpha=0.5)
    axes[1].set_xticks(xq); axes[1].set_xticklabels(lbls)
    axes[1].set_ylabel("LOSO-CV mean AUC (k=t+2)")
    axes[1].set_title("B. AUC by |dφ| quartile (Prediction 1)")
    axes[1].set_ylim(0.45, 0.85)
    axes[1].grid(axis="y", alpha=0.3); axes[1].legend(loc="lower right", fontsize=9)

    plt.suptitle("SIC-prediction experiment (LOSO-CV over 26 sessions, "
                 f"N={m['n_points_total']} points, events={m['n_events_total']})",
                 y=1.02)
    plt.tight_layout()
    p3 = os.path.join(REPORTS, "sic_prediction_leadtime.png")
    plt.savefig(p3, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", p3)


if __name__ == "__main__":
    main()
