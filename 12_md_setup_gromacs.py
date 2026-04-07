#!/usr/bin/env python3
"""
Set up GROMACS simulations for each protein-Nd³⁺ complex.
"""
import subprocess
import os
import yaml

def run_gmx_pdb2gmx(pdb_file, force_field="charmm36-jul2022.ff"):
    # Use pdb2gmx with custom residue for Nd? Simpler: use existing parameter.
    # For metal ions, we need to add manually. Placeholder.
    cmd = f"gmx pdb2gmx -f {pdb_file} -o processed.gro -p topol.top -ff {force_field} -water tip3p"
    subprocess.run(cmd, shell=True, check=True)
    return "processed.gro", "topol.top"

def add_nd_ion(gro_file, top_file, nd_position="0.000 0.000 0.000"):
    # Manually edit files to include Nd³⁺
    # Placeholder: create a new gro line
    with open(gro_file, "r") as f:
        lines = f.readlines()
    # Add Nd at end (simplified)
    with open(gro_file, "a") as f:
        f.write(f"    1ND     ND     1    {nd_position}\n")
    with open(top_file, "a") as f:
        f.write("; Include Nd³⁺ parameters\n#include \"nd.itp\"\n")
    print("Nd³⁺ added manually. Ensure nd.itp is present.")
    return gro_file, top_file

def solvate_and_ions(gro_file, top_file, box_size=1.2):
    subprocess.run(f"gmx editconf -f {gro_file} -o box.gro -c -d {box_size} -bt cubic", shell=True)
    subprocess.run("gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top", shell=True)
    subprocess.run("gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr", shell=True)
    subprocess.run("echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral", shell=True)
    return "solv_ions.gro", "topol.top"

def main():
    config = yaml.safe_load(open("config.yaml"))
    for target in config["md_targets"]:
        os.makedirs(f"outputs/md/{target['name']}", exist_ok=True)
        os.chdir(f"outputs/md/{target['name']}")
        gro, top = run_gmx_pdb2gmx(f"../../../data/pdbs/{target['pdb']}.pdb")
        gro, top = add_nd_ion(gro, top)
        gro, top = solvate_and_ions(gro, top)
        # Copy necessary mdp files
        subprocess.run("cp ../../../mdp_files/*.mdp .", shell=True)
        os.chdir("../../..")
    print("MD setup completed.")

if __name__ == "__main__":
    main()