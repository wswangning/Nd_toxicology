#!/usr/bin/env python3
"""
Construct PPI networks using STRING API for each organ intersection set.
Save as GraphML and edge lists.
"""
import requests
import pandas as pd
import networkx as nx
import time

def query_string(genes, species=9606, score_threshold=700):
    """Query STRING API for interaction network"""
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "\r".join(genes),
        "species": species,
        "required_score": score_threshold,
        "add_nodes": 50
    }
    response = requests.post(url, data=params)
    data = response.json()
    edges = []
    for item in data:
        edges.append((item["preferredName_A"], item["preferredName_B"], item["score"]))
    return edges

def build_graph(edges):
    G = nx.Graph()
    for u, v, w in edges:
        G.add_edge(u, v, weight=w/1000.0)
    return G

def main():
    organs = ["Hepatotoxicity", "Nephrotoxicity", "Neurotoxicity", "Cardiotoxicity"]
    for organ in organs:
        df = pd.read_csv(f"outputs/{organ}_intersection.csv", header=None)
        genes = df[0].tolist()
        if len(genes) < 2:
            print(f"Not enough genes for {organ}")
            continue
        edges = query_string(genes)
        G = build_graph(edges)
        nx.write_graphml(G, f"outputs/{organ}_ppi.graphml")
        # save edge list
        nx.write_edgelist(G, f"outputs/{organ}_ppi_edges.txt", data=["weight"])
        print(f"Saved {organ} PPI: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

if __name__ == "__main__":
    main()