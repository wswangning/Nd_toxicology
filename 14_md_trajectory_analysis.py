#!/usr/bin/env python3
"""
Compute RMSD, RMSF, hydrogen bonds, distances.
"""
import subprocess
import os
import numpy as np

def calc_rmsd(xtc, tpr, output):
    cmd = f"echo 4 4 | gmx rms -s {tpr} -f {xtc} -o {output}_rmsd.xvg -fit rot+trans"
    subprocess.run(cmd, shell=True, check=True)

def calc_rmsf(xtc, tpr, output):
    cmd = f"echo 4 | gmx rmsf -s {tpr} -f {xtc} -o {output}_rmsf.xvg -res"
    subprocess.run(cmd, shell=True, check=True)

def main():
    targets = ["Nrf2", "GPX4", "SLC7A11", "MT1", "MT2", "ACSL4"]
    for target in targets:
        os.chdir(f"outputs/md/{target}")
        calc_rmsd("md.xtc", "md.tpr", target)
        calc_rmsf("md.xtc", "md.tpr", target)
        # Also distance between Nd and key residue
        # Use gmx distance
        os.chdir("../../..")
    print("Analysis completed.")

if __name__ == "__main__":
    main()