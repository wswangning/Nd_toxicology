#!/usr/bin/env python3
"""
Calculate binding free energy using gmx_MMPBSA.
"""
import subprocess
import os

def run_mmpbsa(tpr, xtc, index, nd_resid):
    # Write input file for gmx_MMPBSA
    inp = f"""&general
sys_name="complex"
forcefield="charmm36"
startframe=5000
endframe=10000
interval=10
/
&gb
igb=5
saltcon=0.150
/
&pb
/
"""
    with open("mmpbsa.in", "w") as f:
        f.write(inp)
    cmd = f"gmx_MMPBSA -O -i mmpbsa.in -cs {tpr} -ct {xtc} -cp {tpr} -rg {nd_resid} -lg {nd_resid}"
    subprocess.run(cmd, shell=True, check=True)
    # Parse results
    os.system("cat FINAL_RESULTS_MMPBSA.dat")

def main():
    targets = ["Nrf2", "GPX4", "SLC7A11", "MT1", "MT2", "ACSL4"]
    for target in targets:
        os.chdir(f"outputs/md/{target}")
        run_mmpbsa("md.tpr", "md.xtc", "index.ndx", "ND")
        os.chdir("../../..")

if __name__ == "__main__":
    main()