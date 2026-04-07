#!/usr/bin/env python3
"""
Parse Vina logs, extract binding affinities, generate interaction diagrams.
"""
import os
import pandas as pd
import glob

def parse_vina_log(logfile):
    affinities = []
    with open(logfile) as f:
        for line in f:
            if line.startswith("    1 "):
                parts = line.split()
                affinities.append(float(parts[1]))
    return affinities

def main():
    logs = glob.glob("outputs/docking/*_log.txt")
    results = []
    for log in logs:
        name = os.path.basename(log).replace("_log.txt", "")
        aff = parse_vina_log(log)
        best = aff[0] if aff else None
        results.append({"target": name, "best_affinity_kcal_mol": best})
    df = pd.DataFrame(results)
    df.to_csv("outputs/docking_binding_affinities.csv", index=False)
    print("Affinities saved.")
    # Optionally call PyMOL to generate figures
    # subprocess.run("pymol docking_vis.pml", shell=True)

if __name__ == "__main__":
    main()