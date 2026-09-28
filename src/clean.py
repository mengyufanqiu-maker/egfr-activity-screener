"""Clean and standardise raw IC50 activity data.

Pipeline:
  1. Drop rows without a usable structure (SMILES) or numeric value.
  2. Standardise every measurement to nanomolar (nM).
  3. Drop implausibly low IC50 values (sub-picomolar, almost certainly
     data-entry or unit errors in the source).
  4. Handle censored values: ``<`` / ``<=`` are kept (they are upper bounds on
     IC50, i.e. the compound is *at least* that potent); ``>`` / ``>=`` are
     dropped because they give no useful upper bound for a potency screen.
  5. Deduplicate: when a molecule was tested many times, keep the single best
     (lowest) IC50.
  6. Convert IC50 (nM) to pIC50 = -log10(IC50 in molar).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config

# ChEMBL unit labels -> multiplier to nanomolar.
UNIT_TO_NM = {"pM": 1e-3, "nM": 1.0, "uM": 1e3, "mM": 1e6, "M": 1e9}


def clean_activities(raw: pd.DataFrame) -> pd.DataFrame:
    df = raw.copy()

    if df.empty:
        return pd.DataFrame(columns=[
            "molecule_chembl_id", "canonical_smiles", "standard_relation",
            "standard_value", "standard_units", "ic50_nm", "pIC50"])

    # 1) require a structure and a numeric activity value
    if "canonical_smiles" not in df.columns:
        raise KeyError("activity payload is missing 'canonical_smiles'; check the ChEMBL schema")
    df = df[df["canonical_smiles"].notna()]
    df = df[df["canonical_smiles"].astype(str).str.strip().ne("")]
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df[df["standard_value"].notna()]

    # 2) standardise units to nM
    df["standard_units"] = df["standard_units"].fillna("nM")
    df = df[df["standard_units"].isin(UNIT_TO_NM)]
    df["ic50_nm"] = df["standard_value"] * df["standard_units"].map(UNIT_TO_NM)

    # 3) drop implausibly low IC50 (sub-picomolar)
    df = df[df["ic50_nm"] >= config.MIN_IC50_NM]

    # 4) drop lower-bound-only records
    df["standard_relation"] = df["standard_relation"].fillna("=")
    df = df[~df["standard_relation"].isin([">", ">="])]

    # 5) deduplicate -> keep the best (lowest) IC50 per molecule
    df = (df.sort_values("ic50_nm", ascending=True)
            .drop_duplicates("molecule_chembl_id", keep="first"))

    # 6) pIC50 (M = nM * 1e-9)
    df["pIC50"] = -np.log10(df["ic50_nm"] * 1e-9)

    keep_cols = ["molecule_chembl_id", "canonical_smiles", "standard_relation",
                 "standard_value", "standard_units", "ic50_nm", "pIC50"]
    return df[keep_cols].reset_index(drop=True)
