# EGFR Kinase Activity Screener

A reproducible **drug-discovery data pipeline** that screens the ChEMBL database
for compounds that are both **highly potent** and **drug-like** against the
epidermal growth factor receptor (EGFR) — then **validates** the screen against
approved drugs, **models** activity with a QSAR model, and **compares** EGFR
with its closest kinase relative, HER2.

EGFR is a canonical oncology target: it drives a large share of non-small-cell
lung cancer (NSCLC) and is the target of three generations of approved
tyrosine-kinase inhibitors (gefitinib → erlotinib → osimertinib). That makes it
an ideal, clinically meaningful case study for an end-to-end data-cleaning and
screening workflow.

## Results at a glance

| Stage | Count |
| --- | --- |
| Raw IC50 records fetched (human EGFR, `CHEMBL203`) | 26,600 |
| Unique compounds after cleaning | 11,859 |
| Potent (pIC50 ≥ 7, i.e. IC50 ≤ 100 nM) | 6,936 |
| **Screened hits (potent + drug-like)** | **4,413** |

### Validation against approved drugs

5 of 6 approved EGFR inhibitors are recovered as hits (gefitinib, erlotinib,
afatinib, osimertinib, dacomitinib), and osimertinib is the single most potent
hit overall (IC50 ≈ 2 pM). The one miss — lapatinib — is a "beyond rule-of-five"
drug (MW 581, logP 6.1) that the Lipinski filter correctly excludes, an
expected and explainable exception. See `data/processed/validation_known_drugs.csv`.

### QSAR model

A random forest on ECFP4 fingerprints predicts pIC50 with **R² = 0.72**,
RMSE = 0.74, MAE = 0.51 (n = 11,859 compounds). See `results/plots/05_qsar_parity.png`.

### Selectivity & structural diversity

Among 2,323 compounds tested against both EGFR and HER2: 1,052 dual, 398
EGFR-selective, 108 HER2-selective. The most common hit scaffold is the
**4-anilinoquinazoline** core shared by gefitinib/erlotinib/lapatinib — direct
evidence that the screen recovers the known EGFR pharmacophore.

## What it does

1. **Fetch** — queries the ChEMBL REST API for every IC50 measurement against
   human EGFR, paginating through the full record set (26,600 records).
2. **Clean** — drops unusable records, standardises units to nM, removes
   physically implausible values, handles censored measurements, deduplicates
   (keeping the best IC50), and converts IC50 to pIC50.
3. **Screen** — potency (pIC50 ≥ 7) + Lipinski's Rule of Five (≤ 1 violation,
   judged by ChEMBL's curated flag).
4. **Validate** — checks that the approved EGFR inhibitors are recovered.
5. **QSAR** — predicts pIC50 from ECFP4 fingerprints with a random forest.
6. **Compare** — EGFR vs HER2 selectivity, plus Murcko-scaffold clustering.
7. **Visualise & export** — plots, and hits exported to CSV/JSON/Excel.

## Dependencies

The core screen uses ChEMBL's own curated molecular properties and therefore
needs **no cheminformatics library**. The QSAR and scaffold analyses
additionally use RDKit and scikit-learn. Install everything with
`pip install -r requirements.txt`.

## Project structure

```
egfr-activity-screener/
├── main.py                 # one-command pipeline (all stages)
├── requirements.txt
├── README.md
├── notebooks/
│   └── analysis.ipynb      # the core pipeline as a narrative notebook
├── src/
│   ├── config.py           # tunable parameters and paths
│   ├── fetch.py            # ChEMBL REST API (parallel, paginated)
│   ├── clean.py            # standardise / dedupe / pIC50
│   ├── filter.py           # Lipinski + potency screen
│   ├── validate.py         # validate against approved drugs
│   ├── qsar.py             # ECFP4 + random forest pIC50 model
│   ├── compare.py          # EGFR/HER2 selectivity + scaffolds
│   └── visualize.py        # exploratory plots
├── tests/
│   └── test_pipeline.py    # smoke test (no network required)
├── data/                   # raw + processed CSVs/JSON/XLSX (created on run)
└── results/plots/          # PNG figures (created on run)
```

## Setup & run

Requires Python 3.9+ and internet access. `python main.py` runs the full
pipeline; the first run takes a few minutes (it fetches ~26,600 activity
records, ~15,000 molecule properties, and trains the QSAR model).

```bash
pip install -r requirements.txt
python main.py
```

Each stage can also be run independently, e.g. `python -m src.qsar` (after the
core pipeline has generated `data/processed/`).

## Tests

The cleaning, filtering, and screening logic is covered by a smoke test that
runs on synthetic data (no network needed):

```bash
python tests/test_pipeline.py
```

## Methodology notes

- **Potency threshold.** pIC50 ≥ 7 means IC50 ≤ 100 nM, a common "highly active" bar.
- **Censored values.** `<` / `<=` are upper bounds and retained; `>` / `>=` are dropped.
- **Implausible floor.** IC50 below 1 pM is treated as a source error and dropped.
- **Deduplication.** Only the single best IC50 per molecule is kept.
- **Drug-likeness authority.** Filtering uses ChEMBL's curated `num_ro5_violations`;
  molecules without it are excluded as unassessable.

## Data-quality & cross-check findings

1. **A physically impossible IC50** (~5 × 10⁻⁹ nM, pIC50 17.3) was removed via a
   0.001 nM sanity floor.
2. **Rule-of-Five recomputation vs. ChEMBL's flag** agrees on 99.7% of the final
   hits (14 of 4,413 differ). The residual is traceable to ChEMBL computing the
   rule on the freebase molecular weight (`mw_freebase`) rather than `full_mwt`,
   and to a different logP estimate around the LogP = 5 boundary.

## Limitations & future work

- Screening uses a single activity type (IC50); a production screen would also
  consider Ki/EC50, assay type, and target selectivity.
- The QSAR model is a first pass; a next step is hyperparameter tuning and
  scaffold-split validation to estimate out-of-distribution performance.
- A natural extension is expanding the comparison to the full ErbB family
  (EGFR / HER2 / HER4) to profile selectivity systematically.
