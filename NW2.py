from Bio import pairwise2
import matplotlib.pyplot as plt
import numpy as np

def needleman_wunsch(seq1, seq2):
    #? m = match score of identical chars
    #? s = open and extend gap penalty for both sequences
    #? parameters = (seq1, seq1, match, mismatch, opening a gap, extending a gap)
    #// scores from https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm 
    alignment = pairwise2.align.globalms(seq1, seq2, 2,-1,-5,-1)
    
    return alignment

def main():
    sequences = [
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAACGAAAAAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAATTAAAAAAAAAAAAAAAAAAAAAAA"
    ]

    num_seqs = len(sequences)

    alignment = []

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            alignment.append(needleman_wunsch(sequences[i], sequences[j])[0].seqA)

    align = list(set(alignment))

    max_len = max([len(i) for i in align])

    for i in align:
        if len(i) < max_len:
            kambing = align.index(i)
            align.pop(kambing)
    
    print(align)

main()