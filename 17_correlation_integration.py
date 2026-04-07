#!/usr/bin/env python3
"""
Correlate computational predictions (docking score, ΔG_bind) with experimental outcomes (EC50, cell death).
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def main():
    docking = pd.read_csv("outputs/docking_binding_affinities.csv")
    mmpbsa = pd.read_csv("outputs/mmpbsa_energies.csv")  # created in step 15
    exp = pd.read_csv("outputs/experimental_summary.csv")  # from step 16
    merged = pd.merge(docking, mmpbsa, on="target")
    merged = pd.merge(merged, exp, on="target")
    # Correlation: binding affinity vs EC50
    corr, p = pearsonr(merged["best_affinity_kcal_mol"], merged["EC50_uM"])
    print(f"Correlation affinity-EC50: r={corr:.3f}, p={p:.3f}")
    plt.scatter(merged["best_affinity_kcal_mol"], merged["EC50_uM"])
    plt.xlabel("Docking affinity (kcal/mol)")
    plt.ylabel("EC50 (μM)")
    plt.savefig("outputs/correlation_docking_ec50.png")
    # Also for ΔG_bind vs viability
    corr2, p2 = pearsonr(merged["deltaG_bind"], merged["viability_5uM"])
    print(f"Correlation ΔG_bind-viability: r={corr2:.3f}, p={p2:.3f}")

if __name__ == "__main__":
    main()