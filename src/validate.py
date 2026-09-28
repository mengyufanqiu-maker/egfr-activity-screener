"""Validate the screen against approved EGFR inhibitors.

A good screen should recover the clinically approved EGFR tyrosine-kinase
inhibitors as hits. This module checks each one against the cleaned dataset and
the hit list, and flags the one expected exception: **lapatinib**, a
"beyond rule-of-five" drug that is potent but (correctly, by design) excluded
by the Lipinski filter because it is large (MW > 500) and lipophilic
(logP > 5).
"""
from __future__ import annotations

import pandas as pd

from . import config

# drug name -> ChEMBL id of the freebase form (max_phase = 4, approved)
KNOWN_DRUGS = {
    "gefitinib": "CHEMBL939",
    "erlotinib": "CHEMBL553",
    "afatinib": "CHEMBL1173655",
    "osimertinib": "CHEMBL3353410",
    "lapatinib": "CHEMBL554",
    "dacomitinib": "CHEMBL2105719",
}

COLS = ["drug", "chembl_id", "pIC50", "ic50_nm", "mwt", "alogp",
        "num_ro5_violations", "in_hits"]


def run(cleaned: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    """Return a validation table for the known EGFR inhibitors."""
    hit_ids = set(hits["molecule_chembl_id"])
    rows = []
    for drug, cid in KNOWN_DRUGS.items():
        rec = cleaned[cleaned["molecule_chembl_id"] == cid]
        if rec.empty:
            rows.append({c: None for c in COLS}
                        | {"drug": drug, "chembl_id": cid, "in_hits": False})
            continue
        r = rec.iloc[0]
        rows.append({
            "drug": drug,
            "chembl_id": cid,
            "pIC50": round(float(r["pIC50"]), 2),
            "ic50_nm": round(float(r["ic50_nm"]), 3),
            "mwt": round(float(r["mwt"]), 1),
            "alogp": round(float(r["alogp"]), 2),
            "num_ro5_violations": int(r["num_ro5_violations"]),
            "in_hits": cid in hit_ids,
        })
    return pd.DataFrame(rows, columns=COLS)


if __name__ == "__main__":
    cleaned = pd.read_csv(config.PROCESSED_DIR / "egfr_ic50_cleaned.csv")
    hits = pd.read_csv(config.PROCESSED_DIR / "egfr_hits.csv")
    table = run(cleaned, hits)
    table.to_csv(config.PROCESSED_DIR / "validation_known_drugs.csv", index=False)
    n_hit = int(table["in_hits"].sum())
    print(f"approved EGFR inhibitors recovered as hits: {n_hit}/{len(table)}")
    print(table.to_string(index=False))
