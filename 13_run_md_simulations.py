#!/usr/bin/env python3
"""
Run energy minimization, NVT, NPT, and production MD (100 ns).
"""
import subprocess
import os

def run_gromacs_step(tpr, output, step_name):
    if step_name == "em":
        cmd = f"gmx mdrun -v -deffnm {output}"
    else:
        cmd = f"gmx mdrun -v -deffnm {output} -cpi {output}.cpt"
    subprocess.run(cmd, shell=True, check=True)

def main():
    targets = ["Nrf2", "GPX4", "SLC7A11", "MT1", "MT2", "ACSL4"]  # example
    for target in targets:
        os.chdir(f"outputs/md/{target}")
        # Energy minimization
        subprocess.run("gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr", shell=True)
        run_gromacs_step("em.tpr", "em", "em")
        # NVT
        subprocess.run("gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr", shell=True)
        run_gromacs_step("nvt.tpr", "nvt", "nvt")
        # NPT
        subprocess.run("gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr", shell=True)
        run_gromacs_step("npt.tpr", "npt", "npt")
        # Production (100 ns)
        subprocess.run("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr", shell=True)
        run_gromacs_step("md.tpr", "md", "md")
        os.chdir("../../..")
    print("MD simulations finished.")

if __name__ == "__main__":
    main()