#!/usr/bin/env python3
"""
Create dotplots, barplots for enrichment results.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_dotplot(df, title, outfile):
    df = df.sort_values("pvalue")
    plt.figure(figsize=(8, 5))
    sns.scatterplot(data=df, y="Term", x=-np.log10(df["pvalue"]), size="Genes", sizes=(20, 200))
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

def main():
    organs = ["Hepatotoxicity", "Nephrotoxicity", "Neurotoxicity", "Cardiotoxicity"]
    for organ in organs:
        kegg = pd.read_csv(f"outputs/{organ}_KEGG.csv")
        if not kegg.empty:
            plot_dotplot(kegg, f"{organ} KEGG", f"outputs/{organ}_KEGG_dotplot.png")
    # Integrated network (chord diagram) requires additional libs; simplified
    print("Visualizations saved.")

if __name__ == "__main__":
    import numpy as np
    main()