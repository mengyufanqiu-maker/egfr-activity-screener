"""Central configuration for the EGFR activity screener.

Every tunable value lives here so the pipeline can be re-run against a
different target or with different thresholds without touching the logic.
"""
from __future__ import annotations

from pathlib import Path

# --- ChEMBL API ---------------------------------------------------------------
CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PAGE_SIZE = 1000            # records per request (ChEMBL caps the limit at 1000)
REQUEST_TIMEOUT = 60        # seconds per request
RETRIES = 3                 # retries with exponential backoff
MAX_WORKERS = 8             # concurrent requests for the I/O-bound fetches
MOLECULE_BATCH_SIZE = 100   # molecule ids per property-lookup request

# --- Target -------------------------------------------------------------------
# Human epidermal growth factor receptor (EGFR / erbB1).
TARGET_PREF_NAME = "Epidermal growth factor receptor erbB1"
TARGET_ORGANISM = "Homo sapiens"
TARGET_CHEMBL_ID_FALLBACK = "CHEMBL203"   # well-known human EGFR id

# --- Activity -----------------------------------------------------------------
STANDARD_TYPE = "IC50"      # screen IC50 inhibition values only

# --- Screening thresholds -----------------------------------------------------
PIC50_CUTOFF = 7.0          # pIC50 >= 7  <=>  IC50 <= 100 nM ("high potency")
MIN_IC50_NM = 0.001         # drop sub-picomolar IC50 as implausible (source errors)
LIPINSKI_MAX_VIOLATIONS = 1  # standard practice: at most one Rule-of-Five violation

# --- Paths --------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results"
PLOTS_DIR = RESULTS_DIR / "plots"

for _d in (DATA_DIR, RAW_DIR, PROCESSED_DIR, RESULTS_DIR, PLOTS_DIR):
    _d.mkdir(parents=True, exist_ok=True)
