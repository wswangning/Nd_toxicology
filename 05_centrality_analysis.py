#!/usr/bin/env python3
"""
Topological analysis: degree, betweenness, closeness centrality.
Identify core targets ( >2*median degree, >median betweenness, >median closeness).
"""
import networkx as nx
import pandas as pd
import numpy as np

def analyze_centralities(G):
    deg = dict(G.degree())
    between = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)
    return deg, between, closeness

def filter_core(deg, between, closeness):
    deg_vals = np.array(list(deg.values()))
    between_vals = np.array(list(between.values()))
    closeness_vals = np.array(list(closeness.values()))
    deg_thresh = 2 * np.median(deg_vals)
    between_thresh = np.median(between_vals)
    closeness_thresh = np.median(closeness_vals)
    core = []
    for node in deg:
        if deg[node] >= deg_thresh and between[node] >= between_thresh and closeness[node] >= closeness_thresh:
            core.append(node)
    return core

def main():
    organs = ["Hepatotoxicity", "Nephrotoxicity", "Neurotoxicity", "Cardiotoxicity"]
    for organ in organs:
        G = nx.read_graphml(f"outputs/{organ}_ppi.graphml")
        deg, between, closeness = analyze_centralities(G)
        core = filter_core(deg, between, closeness)
        pd.Series(core).to_csv(f"outputs/{organ}_core_targets.csv", index=False)
        print(f"{organ} core targets: {len(core)}")

if __name__ == "__main__":
    main()