"""Smoke tests for the EGFR activity screener pipeline.

Runs the cleaning -> Lipinski -> screening -> visualisation stages against a
small synthetic dataset that mirrors the ChEMBL response schema, so the logic
can be verified without network access. Runnable directly:

    python tests/test_pipeline.py

(Also discoverable by pytest if it is installed.)
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from src import clean, config, filter, visualize


def _make_raw() -> pd.DataFrame:
    """Synthetic activity records in the ChEMBL /activity.json schema."""
    rows = [
        dict(molecule_chembl_id="CHEMBL1", canonical_smiles="CCO", standard_relation="=", standard_value=50.0, standard_units="nM"),
        dict(molecule_chembl_id="CHEMBL1", canonical_smiles="CCO", standard_relation="=", standard_value=500.0, standard_units="nM"),  # duplicate -> keep best
        dict(molecule_chembl_id="CHEMBL2", canonical_smiles="CCN", standard_relation="<", standard_value=100.0, standard_units="nM"),  # censored, kept
        dict(molecule_chembl_id="CHEMBL3", canonical_smiles="CCC", standard_relation="=", standard_value=0.5, standard_units="uM"),   # 0.5 uM -> 500 nM (not potent)
        dict(molecule_chembl_id="CHEMBL4", canonical_smiles="CCCl", standard_relation=">", standard_value=10.0, standard_units="nM"),  # lower bound -> dropped
        dict(molecule_chembl_id="CHEMBL5", canonical_smiles="CCCCCC", standard_relation="=", standard_value=5.0, standard_units="nM"),  # potent but 4 Ro5 violations
        dict(molecule_chembl_id="CHEMBL6", canonical_smiles="CCCC", standard_relation="=", standard_value=20.0, standard_units="nM"),  # potent, 1 violation (allowed)
        dict(molecule_chembl_id="CHEMBL7", canonical_smiles="CC", standard_relation="=", standard_value=0.000000005012, standard_units="nM"),  # absurd sub-picomolar -> dropped
    ]
    return pd.DataFrame(rows)


def _make_props() -> pd.DataFrame:
    """Synthetic molecule records in the ChEMBL /molecule.json schema."""
    def mp(mwt, alogp, hbd, hba, psa, ro5):
        return dict(full_mwt=mwt, alogp=alogp, hbd=hbd, hba=hba, psa=psa, num_ro5_violations=ro5)

    rows = [
        dict(molecule_chembl_id="CHEMBL1", molecule_properties=mp(400, 3.0, 1, 5, 80, 0)),
        dict(molecule_chembl_id="CHEMBL2", molecule_properties=mp(300, 2.0, 2, 4, 70, 0)),
        dict(molecule_chembl_id="CHEMBL3", molecule_properties=mp(350, 2.5, 1, 4, 75, 0)),
        dict(molecule_chembl_id="CHEMBL5", molecule_properties=mp(600, 6.0, 6, 12, 150, 4)),
        dict(molecule_chembl_id="CHEMBL6", molecule_properties=mp(520, 3.0, 2, 4, 90, 1)),
    ]
    return pd.DataFrame(rows)


def test_clean() -> pd.DataFrame:
    cleaned = clean.clean_activities(_make_raw())
    assert len(cleaned) == 5, f"expected 5 unique compounds, got {len(cleaned)}"
    assert "CHEMBL4" not in set(cleaned["molecule_chembl_id"]), "'>' record should be dropped"
    assert "CHEMBL7" not in set(cleaned["molecule_chembl_id"]), "implausible sub-picomolar IC50 should be dropped"

    c1 = cleaned[cleaned["molecule_chembl_id"] == "CHEMBL1"].iloc[0]
    assert abs(c1["ic50_nm"] - 50.0) < 1e-9, "dedup should keep the best (lowest) IC50"
    assert abs(c1["pIC50"] - 7.301) < 1e-3

    c3 = cleaned[cleaned["molecule_chembl_id"] == "CHEMBL3"].iloc[0]
    assert abs(c3["ic50_nm"] - 500.0) < 1e-9, "0.5 uM should convert to 500 nM"
    print("  clean: OK")
    return cleaned


def test_filter_and_screen() -> tuple[pd.DataFrame, pd.DataFrame]:
    cleaned = clean.clean_activities(_make_raw())
    merged = filter.add_lipinski(cleaned, _make_props())

    mismatch = merged[merged["lipinski_violations"] != merged["num_ro5_violations"]]
    assert mismatch.empty, f"recomputed Lipinski violations disagree with ChEMBL:\n{mismatch}"

    hits = filter.screen(merged)
    assert len(hits) == 3, f"expected 3 hits, got {len(hits)}"
    assert list(hits["molecule_chembl_id"]) == ["CHEMBL6", "CHEMBL1", "CHEMBL2"], "hits should be sorted by potency"
    print("  filter/screen: OK")
    return merged, hits


def test_visualize() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        config.PLOTS_DIR = Path(tmp)
        cleaned = clean.clean_activities(_make_raw())
        merged = filter.add_lipinski(cleaned, _make_props())
        hits = filter.screen(merged)
        visualize.run(merged, hits, raw_count=len(_make_raw()))
        pngs = sorted(config.PLOTS_DIR.glob("*.png"))
        assert len(pngs) == 4, f"expected 4 plots, got {[p.name for p in pngs]}"
    print("  visualize: OK")


def main() -> None:
    test_clean()
    test_filter_and_screen()
    test_visualize()
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
