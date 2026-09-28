"""Fetch EGFR bioactivity data from the ChEMBL REST API.

Three operations are performed:
  1. Resolve the ChEMBL target id for human EGFR (fallback: CHEMBL203).
  2. Page through every IC50 activity record for that target.
  3. Look up ChEMBL's curated molecular properties for the unique molecules.

Requests are issued through a small thread pool because the calls are
independent, I/O-bound HTTP GETs — parallelism cuts the total runtime from
minutes to seconds. No external chemistry library is required: Lipinski
descriptors come straight from ChEMBL's own `molecule_properties`, keeping our
thresholds consistent with the database's own Rule-of-Five flag.
"""
from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import pandas as pd
import requests

from . import config


def _get_json(url: str, params: dict[str, Any] | None = None) -> dict:
    """GET a ChEMBL endpoint and return parsed JSON, with retry/backoff."""
    params = params or {}
    last_err: Exception | None = None
    for attempt in range(config.RETRIES):
        try:
            resp = requests.get(url, params=params, timeout=config.REQUEST_TIMEOUT)
            resp.raise_for_status()
            return resp.json()
        except Exception as e:  # network / HTTP / JSON decoding errors
            last_err = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed to fetch {url} after {config.RETRIES} attempts: {last_err}")


def resolve_target(pref_name: str = config.TARGET_PREF_NAME) -> str:
    """Return the ChEMBL target id for human EGFR.

    Searches by preferred name, filters to *Homo sapiens*, then falls back to
    the well-known human EGFR id (CHEMBL203) if nothing clean matches.
    """
    url = f"{config.CHEMBL_API_BASE}/target/search.json"
    data = _get_json(url, params={"q": pref_name, "limit": 50})
    for t in data.get("targets", []):
        pn = (t.get("pref_name") or "").lower()
        if (t.get("organism") == config.TARGET_ORGANISM
                and "epidermal growth factor" in pn
                and "receptor" in pn):
            return t["target_chembl_id"]
    return config.TARGET_CHEMBL_ID_FALLBACK


def fetch_activities(target_chembl_id: str,
                     standard_type: str = config.STANDARD_TYPE) -> pd.DataFrame:
    """Page through all activity records for the target and return a DataFrame."""
    url = f"{config.CHEMBL_API_BASE}/activity.json"
    base_params = {
        "target_chembl_id": target_chembl_id,
        "standard_type": standard_type,
        "limit": config.PAGE_SIZE,
    }

    # First page also gives us the total record count to paginate against.
    first = _get_json(url, params={**base_params, "offset": 0})
    total = int(first.get("page_meta", {}).get("total_count", 0))
    rows = list(first.get("activities", []))
    offsets = list(range(config.PAGE_SIZE, total, config.PAGE_SIZE))

    def _fetch_page(offset: int) -> list[dict]:
        return _get_json(url, params={**base_params, "offset": offset}).get("activities", [])

    with ThreadPoolExecutor(max_workers=config.MAX_WORKERS) as ex:
        for batch in ex.map(_fetch_page, offsets):
            rows.extend(batch)

    return pd.DataFrame(rows)


def fetch_molecule_properties(molecule_ids: list[str]) -> pd.DataFrame:
    """Fetch curated molecular properties for a list of ChEMBL molecule ids.

    Queries are batched (to keep URLs a reasonable length) and issued in
    parallel.
    """
    url = f"{config.CHEMBL_API_BASE}/molecule.json"
    ids = list(dict.fromkeys(molecule_ids))  # dedupe, preserve order
    chunks = [ids[i:i + config.MOLECULE_BATCH_SIZE]
              for i in range(0, len(ids), config.MOLECULE_BATCH_SIZE)]

    def _fetch_chunk(chunk: list[str]) -> list[dict]:
        data = _get_json(url, params={
            "molecule_chembl_id__in": ",".join(chunk),
            "limit": len(chunk),
        })
        return data.get("molecules", [])

    rows: list[dict] = []
    with ThreadPoolExecutor(max_workers=config.MAX_WORKERS) as ex:
        for batch in ex.map(_fetch_chunk, chunks):
            rows.extend(batch)

    return pd.DataFrame(rows)
