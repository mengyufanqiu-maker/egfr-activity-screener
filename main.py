"""EGFR Activity Screener — end-to-end pipeline.

Run:  python main.py

Flow:
    resolve target -> fetch IC50 activities -> fetch molecule properties ->
    clean -> Lipinski + potency screen -> visualise -> export (CSV/JSON/XLSX).
"""
from __future__ import annotations

from src import config, fetch, clean, filter as screen, visualize


def main() -> None:
    print("Resolving human EGFR target...")
    target_id = fetch.resolve_target()
    print(f"  -> target: {target_id}")

    print("Fetching IC50 activity records (paginated)...")
    raw = fetch.fetch_activities(target_id)
    print(f"  -> {len(raw):,} raw records")
    raw.to_csv(config.RAW_DIR / "egfr_ic50_raw.csv", index=False)

    print("Fetching curated molecular properties...")
    molecule_ids = raw["molecule_chembl_id"].dropna().unique().tolist()
    props = fetch.fetch_molecule_properties(molecule_ids)
    print(f"  -> properties for {len(props):,} molecules")

    print("Cleaning and standardising...")
    cleaned = clean.clean_activities(raw)
    print(f"  -> {len(cleaned):,} unique compounds after cleaning")

    print("Applying Lipinski + potency screen...")
    merged = screen.add_lipinski(cleaned, props)
    hits = screen.screen(merged)
    print(f"  -> {len(hits):,} hits (IC50 <= 100 nM & Lipinski <= 1 violation)")

    print("Visualising...")
    visualize.run(merged, hits, raw_count=len(raw))

    print("Exporting...")
    merged.to_csv(config.PROCESSED_DIR / "egfr_ic50_cleaned.csv", index=False)
    hits.to_csv(config.PROCESSED_DIR / "egfr_hits.csv", index=False)
    hits.to_json(config.PROCESSED_DIR / "egfr_hits.json", orient="records", indent=2)
    hits.to_excel(config.PROCESSED_DIR / "egfr_hits.xlsx", index=False)

    print("\nDone. Outputs written to:")
    print(f"  raw:       {config.RAW_DIR}")
    print(f"  processed: {config.PROCESSED_DIR}")
    print(f"  plots:     {config.PLOTS_DIR}")


if __name__ == "__main__":
    main()
