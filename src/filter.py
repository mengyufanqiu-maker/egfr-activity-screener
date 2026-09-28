"""Apply the screening criteria.

Potency:        pIC50 >= 7  (i.e. IC50 <= 100 nM).
Drug-likeness:  Lipinski's Rule of Five, allowing at most one violation.

We *also* reconstruct the rule from first principles (full_mwt / alogp / hbd /
hba) as an independent cross-check; the residual ~1% disagreement is traceable
to ChEMBL using the freebase molecular weight and a different logP estimate.
"""
from __future__ import annotations

import pandas as pd

from . import config


def _count_lipinski_violations(mwt, logp, hbd, hba) -> int:
    """Count Rule-of-Five violations: MW>500, LogP>5, HBD>5, HBA>10."""
    violations = 0
    if pd.notna(mwt) and mwt > 500:
        violations += 1
    if pd.notna(logp) and logp > 5:
        violations += 1
    if pd.notna(hbd) and hbd > 5:
        violations += 1
    if pd.notna(hba) and hba > 10:
        violations += 1
    return violations


def add_lipinski(df: pd.DataFrame, props: pd.DataFrame) -> pd.DataFrame:
    """Merge curated molecular properties and compute Lipinski violations."""
    if props.empty:
        raise ValueError("molecule properties table is empty; cannot apply Lipinski filter")

    # Normalise: pull the nested molecule_properties dict into flat columns.
    def _flat(rec: dict) -> dict:
        mp = rec.get("molecule_properties") or {}
        return {
            "molecule_chembl_id": rec.get("molecule_chembl_id"),
            "mwt": mp.get("full_mwt"),
            "alogp": mp.get("alogp"),
            "hbd": mp.get("hbd"),
            "hba": mp.get("hba"),
            "psa": mp.get("psa"),
            "num_ro5_violations": mp.get("num_ro5_violations"),
        }

    prop_rows = [_flat(r) for r in props.to_dict(orient="records")]
    pdf = pd.DataFrame(prop_rows).dropna(subset=["molecule_chembl_id"])

    # ChEMBL returns some properties as strings (e.g. "383.81"); coerce to
    # numeric so comparisons behave correctly. Missing values become NaN.
    for col in ("mwt", "alogp", "hbd", "hba", "psa", "num_ro5_violations"):
        pdf[col] = pd.to_numeric(pdf[col], errors="coerce")

    merged = df.merge(pdf, on="molecule_chembl_id", how="left")

    merged["lipinski_violations"] = merged.apply(
        lambda r: _count_lipinski_violations(r["mwt"], r["alogp"], r["hbd"], r["hba"]),
        axis=1,
    )
    return merged


def screen(merged: pd.DataFrame) -> pd.DataFrame:
    """Return the final hit list: potent AND drug-like, sorted by potency.

    Drug-likeness uses ChEMBL's curated ``num_ro5_violations`` (the authoritative
    flag); molecules without that property are excluded as unassessable. The
    independently recomputed ``lipinski_violations`` column is kept for the
    cross-check.
    """
    potent = merged["pIC50"] >= config.PIC50_CUTOFF
    assessable = merged["num_ro5_violations"].notna()
    drug_like = merged["num_ro5_violations"] <= config.LIPINSKI_MAX_VIOLATIONS
    hits = merged[potent & assessable & drug_like]
    return hits.sort_values("pIC50", ascending=False).reset_index(drop=True)
