import numpy as np
from Bio import pairwise2
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import time
import tracemalloc

def smith_waterman_distance(seq1, seq2):
    #? m = match score of identical chars
    #? s = open and extend gap penalty for both sequences
    #? parameters = (seq1, seq1, match, mismatch, opening a gap, extending a gap)
    #// scores from https://gist.github.com/radaniba/11019717 
    alignment = pairwise2.align.localms(seq1, seq2, 2, -1, -1, -1)
    #? gap score
    score = alignment[0].score
    return len(seq1) + len(seq2) - score * 2

def main():
    tracemalloc.start()
    start = time.time()
    #// Example sequences
    sequences = [
        "ACCGTACGGA",
        "ACC-ACAGGA",
        "ACCGTA--GA",
        "AGTACGCA",
        "TATGC",
        "AGTACCGA",
        "AGTACCCA"
    ]

    num_seqs = len(sequences)
    distance_matrix = np.zeros((num_seqs, num_seqs))

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            distance = smith_waterman_distance(sequences[i], sequences[j])
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    #// Perform hierarchical clustering (UPGMA)
    link = hierarchy.linkage(distance_matrix, method='average')
    end = time.time()
    mem = tracemalloc.get_traced_memory()
    print(f"Time = {end-start}")
    print(f"memory = {mem[0]}")

    #// Plot the phylogenetic tree
    #// https://codinginfinite.com/plot-dendrogram-in-python/
    plt.figure(figsize=(8, 6))
    dendrogram = hierarchy.dendrogram(link, labels=sequences, leaf_rotation=-75, leaf_font_size=6)
    plt.xlabel('Sequence Index')
    plt.ylabel('Distance')
    plt.title('Phylogenetic Tree')
    plt.show()

main()