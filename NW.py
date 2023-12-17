from Bio import Entrez
from Bio import SeqIO
from Bio import pairwise2
import matplotlib.pyplot as plt
import tracemalloc
import time
import pylab
from scipy.cluster.hierarchy import linkage, dendrogram
import numpy as np

def needleman_wunsch(seq1, seq2):
    alignment = pairwise2.align.globalms(seq1, seq2, 2,-1,-5,-1)
    score = alignment[0].score
    return len(seq1) + len(seq2) - score * 2

def main():
    tracemalloc.start()
    start = time.time()
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
            distance = needleman_wunsch(sequences[i], sequences[j])
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    link = linkage(distance_matrix, method='average')
    end = time.time()
    mem = tracemalloc.get_traced_memory()
    print(f"Time = {end - start}")
    print(f"Memory = {mem[0]}")

    plt.figure(figsize=(8,6))
    PT = dendrogram(link, labels=sequences, leaf_rotation=-75, leaf_font_size=6)
    plt.xlabel("Sequence Index")
    plt.ylabel("Distance")
    plt.title("Phylogenetic Tree")
    plt.show()

main()