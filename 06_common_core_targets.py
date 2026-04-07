#!/usr/bin/env python3
"""
Find common core targets across at least two organ systems.
"""
import pandas as pd

def main():
    organs = ["Hepatotoxicity", "Nephrotoxicity", "Neurotoxicity", "Cardiotoxicity"]
    core_sets = []
    for organ in organs:
        df = pd.read_csv(f"outputs/{organ}_core_targets.csv", header=None)
        core_sets.append(set(df[0].tolist()))
    
    # Count occurrences across organs
    all_genes = set.union(*core_sets)
    counts = {gene: sum(1 for s in core_sets if gene in s) for gene in all_genes}
    common = [gene for gene, cnt in counts.items() if cnt >= 2]
    common_df = pd.DataFrame({"gene": common, "organ_count": [counts[g] for g in common]})
    common_df.to_csv("outputs/common_core_targets.csv", index=False)
    print("Common core targets:", common)

if __name__ == "__main__":
    main()