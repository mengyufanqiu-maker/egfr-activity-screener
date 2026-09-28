"""Build a QSAR model to predict EGFR pIC50 from molecular structure.

Features: ECFP4 (Morgan, radius 2) fingerprints computed with RDKit.
Model:    random-forest regressor (scikit-learn).
Outputs:  R2 / RMSE / MAE on a held-out test set, plus a parity plot.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split

from . import config

N_BITS = 2048          # ECFP4 folded to 2048 bits
RADIUS = 2
TEST_SIZE = 0.2
N_ESTIMATORS = 300
RANDOM_STATE = 42

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=N_BITS)


def run(cleaned: pd.DataFrame) -> dict:
    """Train a random forest on ECFP4 fingerprints and return metrics."""
    X_rows, y_rows = [], []
    n_failed = 0
    for smi, pic in zip(cleaned["canonical_smiles"], cleaned["pIC50"]):
        if not isinstance(smi, str) or not smi.strip():
            n_failed += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_failed += 1
            continue
        fp = _GEN.GetFingerprint(mol)
        X_rows.append(list(fp))
        y_rows.append(float(pic))

    X = np.asarray(X_rows, dtype=np.float32)
    y = np.asarray(y_rows, dtype=np.float64)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=RANDOM_STATE)

    model = RandomForestRegressor(
        n_estimators=N_ESTIMATORS, n_jobs=-1, random_state=RANDOM_STATE)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    return {
        "n_compounds": len(y),
        "n_failed_smiles": n_failed,
        "r2": float(r2_score(y_test, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_test, y_pred))),
        "mae": float(mean_absolute_error(y_test, y_pred)),
        "y_test": y_test,
        "y_pred": y_pred,
        "model": model,
    }


def plot_parity(y_test: np.ndarray, y_pred: np.ndarray, r2: float, rmse: float) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_test, y_pred, s=6, alpha=0.25, color="#2c7fb8")
    lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
    ax.plot(lims, lims, "r--", lw=1.5, label="y = x")
    ax.set_xlabel("Measured pIC50")
    ax.set_ylabel("Predicted pIC50")
    ax.set_title(f"QSAR parity plot (R² = {r2:.2f}, RMSE = {rmse:.2f})")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(config.PLOTS_DIR / "05_qsar_parity.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    cleaned = pd.read_csv(config.PROCESSED_DIR / "egfr_ic50_cleaned.csv")
    res = run(cleaned)
    print(f"compounds modelled: {res['n_compounds']}  (skipped: {res['n_failed_smiles']})")
    print(f"R²   = {res['r2']:.3f}")
    print(f"RMSE = {res['rmse']:.3f}  (in pIC50 units)")
    print(f"MAE  = {res['mae']:.3f}")
    plot_parity(res["y_test"], res["y_pred"], res["r2"], res["rmse"])
    print("parity plot ->", config.PLOTS_DIR / "05_qsar_parity.png")
