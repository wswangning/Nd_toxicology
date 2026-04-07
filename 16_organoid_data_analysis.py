#!/usr/bin/env python3
"""
Process in vitro data: viability, GSH/GSSG, MDA, qPCR, Western blot, flow cytometry.
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

def load_data():
    # Assume CSV files are in data/experimental/
    viability = pd.read_csv("data/experimental/viability.csv")
    gsh = pd.read_csv("data/experimental/gsh_gssg.csv")
    mda = pd.read_csv("data/experimental/mda.csv")
    qpcr = pd.read_csv("data/experimental/qpcr.csv")
    wb = pd.read_csv("data/experimental/western.csv")
    return viability, gsh, mda, qpcr, wb

def stats_analysis(df, control_col, treat_cols):
    results = {}
    for col in treat_cols:
        t, p = ttest_ind(df[control_col], df[col])
        results[col] = p
    return results

def plot_bar(data, labels, title, outfile):
    plt.figure()
    plt.bar(labels, data)
    plt.title(title)
    plt.savefig(outfile, dpi=300)
    plt.close()

def main():
    viability, gsh, mda, qpcr, wb = load_data()
    # Example: plot viability
    doses = [0, 1, 2.5, 5]
    mean_viab = viability.mean()
    plot_bar(mean_viab, doses, "Viability (%)", "outputs/viability.png")
    # Save stats
    stats = stats_analysis(viability, "0uM", ["1uM","2.5uM","5uM"])
    pd.Series(stats).to_csv("outputs/viability_stats.csv")
    print("Organoid data processed.")

if __name__ == "__main__":
    main()