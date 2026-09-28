"""Exploratory plots for the EGFR activity screener.

Produces a small set of publication-style PNGs in results/plots/:
  01_pic50_distribution.png    — pIC50 histogram with the potency cutoff
  02_mwt_vs_pic50.png          — molecular weight vs potency, hits highlighted
  03_lipinski_violations.png   — Rule-of-Five violation counts
  04_screening_funnel.png      — how many compounds survive each stage
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from . import config

sns.set_theme(style="whitegrid", context="notebook")
HIT_COLOR = "#c0392b"
PASS_COLOR = "#2c7fb8"


def run(merged: pd.DataFrame, hits: pd.DataFrame, raw_count: int) -> None:
    hit_ids = set(hits["molecule_chembl_id"])
    merged = merged.assign(is_hit=merged["molecule_chembl_id"].isin(hit_ids))

    # 1) pIC50 distribution ------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.histplot(merged["pIC50"], bins=50, kde=True, color=PASS_COLOR, ax=ax)
    ax.axvline(config.PIC50_CUTOFF, color=HIT_COLOR, ls="--", lw=2,
               label=f"potency cutoff (pIC50 = {config.PIC50_CUTOFF})")
    ax.set_xlabel("pIC50  (−log10 IC50, M)")
    ax.set_ylabel("Number of compounds")
    ax.set_title("Distribution of EGFR potency (IC50)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(config.PLOTS_DIR / "01_pic50_distribution.png", dpi=150)
    plt.close(fig)

    # 2) molecular weight vs potency --------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    non_hits = merged[~merged["is_hit"]]
    ax.scatter(non_hits["mwt"], non_hits["pIC50"], s=12, alpha=0.35,
               color=PASS_COLOR, label="filtered out")
    ax.scatter(hits["mwt"], hits["pIC50"], s=16, alpha=0.8,
               color=HIT_COLOR, label="screened hits")
    ax.axvline(500, color="gray", ls=":", lw=1.5, label="Lipinski MW ≤ 500")
    ax.axhline(config.PIC50_CUTOFF, color=HIT_COLOR, ls="--", lw=1.5)
    ax.set_xlabel("Molecular weight (Da)")
    ax.set_ylabel("pIC50")
    ax.set_title("Molecular weight vs. potency")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(config.PLOTS_DIR / "02_mwt_vs_pic50.png", dpi=150)
    plt.close(fig)

    # 3) Lipinski violation counts ----------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    max_v = int(merged["lipinski_violations"].max()) if merged["lipinski_violations"].notna().any() else 4
    bins = np.arange(-0.5, max(max_v, 4) + 1.5, 1.0)
    counts, _ = np.histogram(merged["lipinski_violations"].fillna(0), bins=bins)
    bar_colors = [PASS_COLOR if i <= config.LIPINSKI_MAX_VIOLATIONS else "#bdc3c7"
                  for i in range(len(bins) - 1)]
    ax.bar(bins[:-1] + 0.5, counts, color=bar_colors)
    ax.axvspan(-0.5, config.LIPINSKI_MAX_VIOLATIONS + 0.5, color=PASS_COLOR, alpha=0.08)
    ax.set_xlabel("Lipinski Rule-of-Five violations")
    ax.set_ylabel("Number of compounds")
    ax.set_title("Drug-likeness (Lipinski) violations")
    ax.set_xticks(range(len(bins) - 1))
    fig.tight_layout()
    fig.savefig(config.PLOTS_DIR / "03_lipinski_violations.png", dpi=150)
    plt.close(fig)

    # 4) screening funnel --------------------------------------------------
    n_cleaned = len(merged)
    n_potent = int((merged["pIC50"] >= config.PIC50_CUTOFF).sum())
    n_hits = len(hits)
    stages = ["Raw records", "Cleaned\n(unique compounds)", "Potent\n(pIC50 ≥ 7)", "Hits\n(potent + drug-like)"]
    values = [raw_count, n_cleaned, n_potent, n_hits]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(stages, values, color=[PASS_COLOR, PASS_COLOR, PASS_COLOR, HIT_COLOR])
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width() / 2, v + max(values) * 0.02, f"{v:,}",
                ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("Count (log scale)")
    ax.set_yscale("log")
    ax.set_title("Screening funnel")
    fig.tight_layout()
    fig.savefig(config.PLOTS_DIR / "04_screening_funnel.png", dpi=150)
    plt.close(fig)
