# import numpy as np
# from Bio import Phylo
# import io

# sequences = [
#     "AGTACGCA",
#     "TATGC",
#     "AGTACCGA",
#     "AGTACCCA"
# ]

# # Smith-Waterman Algorithm for local sequence alignment
# def smith_waterman(seq1, seq2, match_score=3, mismatch_score=-1, gap_penalty=-2):
#     # Initialization
#     m, n = len(seq1), len(seq2)
#     score_matrix = np.zeros((m + 1, n + 1))
#     traceback_matrix = np.zeros((m + 1, n + 1))

#     # Fill the matrices
#     for i in range(1, m + 1):
#         for j in range(1, n + 1):
#             match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
#             delete = score_matrix[i - 1][j] + gap_penalty
#             insert = score_matrix[i][j - 1] + gap_penalty
#             score_matrix[i][j] = max(0, match, delete, insert)
#             if score_matrix[i][j] == 0:
#                 traceback_matrix[i][j] = 0  # 0 indicates the end of the local alignment
#             elif score_matrix[i][j] == match:
#                 traceback_matrix[i][j] = 1  # 1 indicates a match
#             elif score_matrix[i][j] == delete:
#                 traceback_matrix[i][j] = 2  # 2 indicates an insertion
#             else:
#                 traceback_matrix[i][j] = 3  # 3 indicates a deletion

#     # Traceback to find the alignment with the highest score
#     max_score = np.max(score_matrix)
#     max_pos = np.where(score_matrix == max_score)
#     i, j = max_pos[0][0], max_pos[1][0]

#     aligned_seq1, aligned_seq2 = '', ''
#     while traceback_matrix[i][j] != 0:
#         if traceback_matrix[i][j] == 1:
#             aligned_seq1 = seq1[i - 1] + aligned_seq1
#             aligned_seq2 = seq2[j - 1] + aligned_seq2
#             i -= 1
#             j -= 1
#         elif traceback_matrix[i][j] == 2:
#             aligned_seq1 = seq1[i - 1] + aligned_seq1
#             aligned_seq2 = '-' + aligned_seq2
#             i -= 1
#         elif traceback_matrix[i][j] == 3:
#             aligned_seq1 = '-' + aligned_seq1
#             aligned_seq2 = seq2[j - 1] + aligned_seq2
#             j -= 1

#     return max_score / max(m, n), aligned_seq1, aligned_seq2  # Returning alignment score and aligned sequences


# # Calculate distances based on alignment scores using Smith-Waterman
# def calculate_distances(sequences):
#     distances = np.zeros((len(sequences), len(sequences)))
#     for i in range(len(sequences)):
#         for j in range(i + 1, len(sequences)):
#             alignment_score, _, _ = smith_waterman(sequences[i], sequences[j])
#             distances[i, j] = alignment_score
#             distances[j, i] = alignment_score
#     return distances


# # UPGMA algorithm
# def upgma(dist_matrix, names):
#     num_seqs = len(dist_matrix)
#     cluster_names = names.copy()
#     tree = {}

#     while len(cluster_names) > 1:
#         min_index = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
#         min_dist = dist_matrix[min_index]

#         # Find the sequences corresponding to the minimum distance
#         seq_idx1, seq_idx2 = min_index

#         # Create a new node
#         new_node = f'({cluster_names[seq_idx1]}:{min_dist/2},{cluster_names[seq_idx2]}:{min_dist/2})'

#         # Update the distance matrix and cluster names
#         new_dist_row = (dist_matrix[seq_idx1] + dist_matrix[seq_idx2]) / 2
#         new_dist_row = np.reshape(new_dist_row, (len(new_dist_row), 1))
#         new_dist_col = np.append(new_dist_row, [[0]], axis=0)  # Add a placeholder for the new node's distance to itself

#         dist_matrix = np.delete(dist_matrix, [seq_idx1, seq_idx2], axis=0)
#         dist_matrix = np.delete(dist_matrix, [seq_idx1, seq_idx2], axis=1)

#         dist_matrix = np.insert(dist_matrix, len(dist_matrix), new_dist_col, axis=1)
#         new_dist_row = np.append(new_dist_row, [[0]], axis=0)  # Add a placeholder for the new node's distance to itself
#         dist_matrix = np.insert(dist_matrix, len(dist_matrix), new_dist_row, axis=0)

#         cluster_names.pop(max(seq_idx1, seq_idx2))  # Remove the larger index first
#         cluster_names.pop(min(seq_idx1, seq_idx2))  # Remove the smaller index next

#         # Store the new node in the tree
#         tree[new_node] = [cluster_names[idx] for idx in range(len(cluster_names))]

#         cluster_names.append(new_node)

#     return tree

# # Calculate distances and construct the tree using UPGMA
# dist_matrix = calculate_distances(sequences)
# tree = upgma(dist_matrix, [f"Seq_{i+1}" for i in range(len(sequences))])

# def draw_phylogenetic_tree(tree):
#     # Convert the tree dictionary into Newick format
#     newick_tree = list(tree.keys())[0]  # Assuming the tree has a single root node
#     newick_tree += ";"

#     # Parse the Newick tree string
#     parsed_tree = Phylo.read(io.StringIO(newick_tree), "newick")

#     # Draw the phylogenetic tree
#     Phylo.draw(parsed_tree)


# # draw_phylogenetic_tree(tree)

# print("Phylogenetic Tree:")
# for node in tree:
#     print(node, ":", tree[node])

import numpy as np
from Bio import pairwise2
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

def smith_waterman_distance(seq1, seq2):
    alignment = pairwise2.align.localms(seq1, seq2, 2, -1, -1, -1)
    score = alignment[0].score
    return len(seq1) + len(seq2) - score * 2

# Example sequences
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

# Perform hierarchical clustering (UPGMA)
linkage = hierarchy.linkage(distance_matrix, method='average')

# Plot the dendrogram
plt.figure(figsize=(8, 6))
dendrogram = hierarchy.dendrogram(linkage, labels=sequences, leaf_rotation=-90)
plt.xlabel('Sequence Index')
plt.ylabel('Distance')
plt.title('Phylogenetic Tree')
plt.show()

