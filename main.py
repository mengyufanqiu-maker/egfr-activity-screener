"""EGFR Activity Screener — end-to-end pipeline.

Run:  python main.py

Stages:
  1. resolve target          — human EGFR (CHEMBL203)
  2. fetch activities        — all IC50 records (paginated, parallel)
  3. fetch properties        — ChEMBL curated molecular properties
  4. clean + screen          — pIC50 + Lipinski (<= 1 violation)
  5. validate                — recover the approved EGFR inhibitors
  6. qsar                    — predict pIC50 from ECFP4 (random forest)
  7. compare                 — EGFR vs HER2 selectivity + Murcko scaffolds
  8. visualise + export      — plots, CSV / JSON / XLSX
"""
from __future__ import annotations

from src import config, fetch, clean, filter as screen, validate, visualize
from src import qsar, compare


def main() -> None:
    print("== 1/8 Resolving human EGFR target ==")
    target_id = fetch.resolve_target()
    print(f"  target: {target_id}")

    print("== 2/8 Fetching IC50 activity records (paginated) ==")
    raw = fetch.fetch_activities(target_id)
    print(f"  {len(raw):,} raw records")
    raw.to_csv(config.RAW_DIR / "egfr_ic50_raw.csv", index=False)

    print("== 3/8 Fetching curated molecular properties ==")
    molecule_ids = raw["molecule_chembl_id"].dropna().unique().tolist()
    props = fetch.fetch_molecule_properties(molecule_ids)
    print(f"  properties for {len(props):,} molecules")

    print("== 4/8 Cleaning and screening ==")
    cleaned = clean.clean_activities(raw)
    merged = screen.add_lipinski(cleaned, props)
    hits = screen.screen(merged)
    print(f"  {len(cleaned):,} unique compounds -> {len(hits):,} hits")
    # save core outputs early so a partial run still leaves usable data
    merged.to_csv(config.PROCESSED_DIR / "egfr_ic50_cleaned.csv", index=False)
    hits.to_csv(config.PROCESSED_DIR / "egfr_hits.csv", index=False)
    hits.to_json(config.PROCESSED_DIR / "egfr_hits.json", orient="records", indent=2)
    hits.to_excel(config.PROCESSED_DIR / "egfr_hits.xlsx", index=False)

    print("== 5/8 Validating against approved EGFR inhibitors ==")
    validation = validate.run(cleaned, hits)
    validation.to_csv(config.PROCESSED_DIR / "validation_known_drugs.csv", index=False)
    print(f"  recovered {int(validation['in_hits'].sum())}/{len(validation)} approved drugs")

    print("== 6/8 QSAR model (ECFP4 + random forest; ~1 min) ==")
    qsar_res = qsar.run(cleaned)
    qsar.plot_parity(qsar_res["y_test"], qsar_res["y_pred"],
                     qsar_res["r2"], qsar_res["rmse"])
    print(f"  R2 = {qsar_res['r2']:.3f}  RMSE = {qsar_res['rmse']:.3f}  MAE = {qsar_res['mae']:.3f}")

    print("== 7/8 EGFR/HER2 selectivity + scaffold clustering ==")
    cmp = compare.run(cleaned, hits)
    print(f"  HER2 compounds: {cmp['her2_n']:,}")
    print(f"  selectivity: {cmp['selectivity']['class'].value_counts().to_dict()}")
    print(f"  top scaffold (n={cmp['scaffolds'].iloc[0]}): {cmp['scaffolds'].index[0]}")

    print("== 8/8 Visualising ==")
    visualize.run(merged, hits, raw_count=len(raw))

    print("\nDone. Outputs written to data/ and results/plots/.")


if __name__ == "__main__":
    main()
