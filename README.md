# EGFR Kinase Activity Screener

A small, reproducible **drug-discovery data pipeline** that screens the ChEMBL
database for compounds that are both **highly potent** and **drug-like**
against the epidermal growth factor receptor (EGFR).

EGFR is a canonical oncology target: it drives a large share of non-small-cell
lung cancer (NSCLC) and is the target of three generations of approved
tyrosine-kinase inhibitors (gefitinib → erlotinib → osimertinib). That makes it
an ideal, clinically meaningful case study for an end-to-end data-cleaning and
screening workflow.

## Results at a glance

| Stage | Count |
| --- | --- |
| Raw IC50 records fetched (human EGFR, `CHEMBL203`) | 26,600 |
| Unique compounds after cleaning | 11,878 |
| Potent (pIC50 ≥ 7, i.e. IC50 ≤ 100 nM) | 6,945 |
| **Screened hits (potent + drug-like)** | **4,413** |

The most potent hit inhibits EGFR with an IC50 of ~2 pM. The full hit list is
exported to CSV, JSON, and Excel under `data/processed/`.

## What it does

1. **Fetch** — queries the ChEMBL REST API for every `IC50` measurement against
   human EGFR (`CHEMBL203`), paginating through the full record set (26,600
   records).
2. **Clean** — drops unusable records, standardises units to nanomolar, removes
   physically implausible values, handles censored measurements, deduplicates
   molecules (keeping the single best IC50), and converts IC50 to
   `pIC50 = -log10(IC50 in M)`.
3. **Filter** — applies two screening criteria:
   - *Potency*: `pIC50 >= 7` (IC50 ≤ 100 nM)
   - *Drug-likeness*: Lipinski's Rule of Five with at most **1** violation,
     judged by ChEMBL's curated `num_ro5_violations` flag.
4. **Visualise** — four exploratory plots (potency distribution,
   molecular-weight vs potency, Lipinski violation counts, screening funnel).
5. **Export** — writes the hit list to CSV, JSON, and Excel.

## Why ChEMBL's properties instead of RDKit?

Rather than re-deriving physicochemical descriptors with a third-party
cheminformatics library, the pipeline uses ChEMBL's own curated
`molecule_properties` (molecular weight, AlogP, H-bond donors/acceptors, polar
surface area, Rule-of-Five violations). This keeps the thresholds consistent
with the database's own standardisation and removes a heavy dependency — while
still demonstrating understanding by recomputing the Rule of Five
independently as a cross-check.

## Project structure

```
egfr-activity-screener/
├── main.py                 # one-command pipeline
├── requirements.txt        # Python dependencies
├── README.md
├── notebooks/
│   └── analysis.ipynb      # the analysis as a narrative notebook
├── src/
│   ├── config.py           # all tunable parameters and paths
│   ├── fetch.py            # ChEMBL REST API access (parallel, paginated)
│   ├── clean.py            # standardise / dedupe / pIC50
│   ├── filter.py           # Lipinski + potency screen
│   └── visualize.py        # exploratory plots
├── tests/
│   └── test_pipeline.py    # smoke test (no network required)
├── data/                   # raw + processed CSVs/JSON/XLSX (created on run)
└── results/plots/          # PNG figures (created on run)
```

## Setup & run

Requires Python 3.9+ and an internet connection (data is fetched live from
ChEMBL). The first run takes a minute or two — the pipeline fetches ~26,600
activity records and looks up properties for ~15,000 molecules, using a small
thread pool to keep it fast.

```bash
# (optional but recommended) create a virtual environment
python -m venv .venv
# Windows:  .venv\Scripts\activate      macOS/Linux:  source .venv/bin/activate

pip install -r requirements.txt
python main.py
```

The pipeline prints the number of raw records, cleaned compounds, and final
hits, then writes:

| Path | Contents |
| --- | --- |
| `data/raw/egfr_ic50_raw.csv` | raw activity records as returned by ChEMBL |
| `data/processed/egfr_ic50_cleaned.csv` | cleaned, deduplicated compounds + properties |
| `data/processed/egfr_hits.csv` / `.json` / `.xlsx` | the screened hit list |
| `results/plots/*.png` | four exploratory figures |

The notebook (`notebooks/analysis.ipynb`) runs the same pipeline live and walks
through the results step by step.

## Tests

The cleaning, filtering, and screening logic is covered by a smoke test that
runs on synthetic data (no network needed):

```bash
python tests/test_pipeline.py
```

## Methodology notes

- **Potency threshold.** `pIC50 >= 7` means the compound inhibits EGFR with an
  IC50 of 100 nM or better — a commonly used bar for "highly active".
- **Censored values.** A `<` or `<=` result is an upper bound (the compound is
  *at least* that potent), so those records are retained; the value is treated
  conservatively. `>` / `>=` give only a lower bound on IC50 and are dropped.
- **Implausible floor.** IC50 values below 1 pM are treated as source errors
  and dropped (see below).
- **Deduplication.** The same molecule is often tested many times across
  assays and publications; only the single best (lowest) IC50 is kept.
- **Drug-likeness authority.** The filter uses ChEMBL's curated
  `num_ro5_violations`; molecules without it are excluded as unassessable.

## Data-quality & cross-check findings

Two data-quality issues surfaced during the analysis and were handled
explicitly:

1. **A physically impossible IC50.** One record reported an IC50 of
   ~5 × 10⁻⁹ nM (pIC50 17.3) — far below the concentration of a single
   molecule. It was removed via a 0.001 nM (1 pM) sanity floor.
2. **Rule-of-Five recomputation vs. ChEMBL's flag.** Independently recomputing
   Lipinski from `full_mwt` / `alogp` / `hbd` / `hba` agrees with ChEMBL's
   curated flag on 99.7% of the final hits (14 of 4,413 differ). The residual
   disagreement is traceable to ChEMBL computing the rule on the *freebase*
   molecular weight (`mw_freebase`) rather than `full_mwt` (which includes
   counterions), and to a slightly different logP estimate around the
   LogP = 5 boundary.

These checks are exactly the kind of "boundary judgement" that separates a
domain-aware screen from a mechanical one.

## Limitations & future work

- Screening is based on a single activity type (IC50); a production screen
  would also consider Ki/EC50, assay type, and target selectivity.
- Using ChEMBL's curated properties avoids a local RDKit install; a natural
  extension is to recompute descriptors locally and compare.
- A promising next step is a QSAR model that predicts `pIC50` from molecular
  fingerprints, turning this screening pipeline into a predictive tool.
