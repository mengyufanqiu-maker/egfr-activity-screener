"""Compare EGFR against its closest kinase relative, HER2 (ERBB2).

Two analyses:
  1. Selectivity: classify compounds tested against both targets as
     EGFR-selective, HER2-selective, dual, or inactive.
  2. Structural diversity: cluster the EGFR hits by Murcko scaffold.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from . import config, fetch, clean

HER2_CHEMBL_ID = "CHEMBL1824"


def run_selectivity(egfr: pd.DataFrame, her2: pd.DataFrame,
                    cutoff: float = config.PIC50_CUTOFF) -> pd.DataFrame:
    """Classify compounds tested against both EGFR and HER2 by potency."""
    e = egfr[["molecule_chembl_id", "pIC50"]].rename(columns={"pIC50": "pIC50_egfr"})
    h = her2[["molecule_chembl_id", "pIC50"]].rename(columns={"pIC50": "pIC50_her2"})
    m = e.merge(h, on="molecule_chembl_id", how="inner")
    m["delta_pIC50"] = m["pIC50_egfr"] - m["pIC50_her2"]
    conds = [
        (m["pIC50_egfr"] >= cutoff) & (m["pIC50_her2"] >= cutoff),
        (m["pIC50_egfr"] >= cutoff) & (m["pIC50_her2"] < cutoff),
        (m["pIC50_egfr"] < cutoff) & (m["pIC50_her2"] >= cutoff),
    ]
    m["class"] = np.select(conds, ["dual", "egfr_selective", "her2_selective"],
                           default="inactive")
    return m


def run_scaffolds(hits: pd.DataFrame, top_n: int = 10) -> pd.Series:
    """Count Murcko scaffolds among the hits, returning the most common."""
    counts: dict[str, int] = {}
    for smi in hits["canonical_smiles"]:
        if not isinstance(smi, str):
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        scaff = MurckoScaffold.GetScaffoldForMol(mol)
        if scaff is None:
            continue
        s = Chem.MolToSmiles(scaff)
        counts[s] = counts.get(s, 0) + 1
    return pd.Series(counts, name="n_hits").sort_values(ascending=False).head(top_n)


def _plot_selectivity(m: pd.DataFrame, cutoff: float) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    colors = {"dual": "#c0392b", "egfr_selective": "#2c7fb8",
              "her2_selective": "#27ae60", "inactive": "#bdc3c7"}
    for cls, grp in m.groupby("class"):
        ax.scatter(grp["pIC50_her2"], grp["pIC50_egfr"], s=5, alpha=0.25,
                   color=colors[cls], label=cls)
    ax.plot([0, 14], [0, 14], "k--", lw=1)
    ax.axhline(cutoff, color="gray", lw=1, ls=":")
    ax.axvline(cutoff, color="gray", lw=1, ls=":")
    ax.set_xlabel("pIC50 (HER2)")
    ax.set_ylabel("pIC50 (EGFR)")
    ax.set_title("EGFR vs HER2 potency")
    ax.legend(loc="upper left", markerscale=3)
    fig.tight_layout()
    fig.savefig(config.PLOTS_DIR / "06_egfr_her2_selectivity.png", dpi=150)
    plt.close(fig)


def run(egfr_cleaned: pd.DataFrame, hits: pd.DataFrame) -> dict:
    """Fetch HER2, run the selectivity + scaffold analyses, return results."""
    her2_raw = fetch.fetch_activities(HER2_CHEMBL_ID)
    her2 = clean.clean_activities(her2_raw)
    m = run_selectivity(egfr_cleaned, her2)
    m.to_csv(config.PROCESSED_DIR / "egfr_her2_selectivity.csv", index=False)
    _plot_selectivity(m, config.PIC50_CUTOFF)
    scaffolds = run_scaffolds(hits)
    scaffolds.to_csv(config.PROCESSED_DIR / "top_scaffolds.csv", header=["n_hits"])
    return {"her2_n": len(her2), "selectivity": m, "scaffolds": scaffolds}


if __name__ == "__main__":
    egfr = pd.read_csv(config.PROCESSED_DIR / "egfr_ic50_cleaned.csv")
    hits = pd.read_csv(config.PROCESSED_DIR / "egfr_hits.csv")
    res = run(egfr, hits)
    print(f"HER2 compounds: {res['her2_n']:,}")
    print("Selectivity:", res["selectivity"]["class"].value_counts().to_dict())
    print("Top scaffolds:")
    for smi, n in res["scaffolds"].items():
        print(f"  {n:5d}  {smi}")
