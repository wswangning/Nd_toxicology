# Neodymium Nitrate Hepatotoxicity – Integrated Computational Framework

This repository contains all custom analysis scripts for the paper:  
**"Mechanistic Insights into Neodymium Nitrate-Induced Hepatotoxicity: An Integrated Framework Combining Network Toxicology, Molecular Dynamics, and 3D Organoid Validation"**  
Wang et al., 2026.

## Overview

We provide Python scripts to reproduce:
- Target screening and PPI network analysis (STRING, Cytoscape)
- GO/KEGG enrichment (DAVID/clusterProfiler)
- Molecular docking (AutoDock Vina)
- Molecular dynamics simulations (GROMACS + gmx_MMPBSA)
- In vitro organoid data analysis (viability, qPCR, Western blot, etc.)
- Correlation between computational predictions and experiments

## Requirements

- Python 3.9+
- Dependencies listed in `requirements.txt` or `environment.yml`
- External software: AutoDock Vina, GROMACS (2022+), gmx_MMPBSA, PyMOL (optional)

## Usage

1. Clone this repository.
2. Set up environment: `conda env create -f environment.yml`
3. Prepare input data in `data/` (see example).
4. Run scripts in numerical order:
   ```bash
   python 01_target_screening.py
   python 02_organ_targets_fetch.py
   ...