"""EGFR activity screener — a small drug-discovery data pipeline.

Modules:
    config      — tunable parameters and paths
    fetch       — ChEMBL REST API access (target, activities, molecule properties)
    clean       — standardise units, handle censored values, deduplicate, pIC50
    filter      — Lipinski's Rule of Five + potency screening
    visualize   — exploratory plots
"""
