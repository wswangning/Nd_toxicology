#!/usr/bin/env python3
"""
Calculate intersections between Nd targets and each organ's disease targets.
Generate Venn diagrams.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os

def load_targets(organ):
    return set(pd.read_csv(f"outputs/{organ}_targets.csv", header=None)[0].tolist())

def main():
    nd_targets = set(pd.read_csv("outputs/nd_targets_raw.csv", header=None)[0].tolist())
    organs = ["Hepatotoxicity", "Nephrotoxicity", "Neurotoxicity", "Cardiotoxicity"]
    intersections = {}
    for organ in organs:
        organ_set = load_targets(organ)
        inter = nd_targets.intersection(organ_set)
        intersections[organ] = inter
        print(f"{organ}: {len(inter)} intersecting targets")
        pd.Series(list(inter)).to_csv(f"outputs/{organ}_intersection.csv", index=False)
        # Venn diagram
        plt.figure()
        venn2([nd_targets, organ_set], set_labels=('Nd targets', organ))
        plt.title(f"{organ} intersection")
        plt.savefig(f"outputs/venn_{organ}.png", dpi=300)
        plt.close()

if __name__ == "__main__":
    main()