import numpy as np
import matplotlib.pyplot as plt

def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-1):
    len_seq1 = len(seq1) + 1
    len_seq2 = len(seq2) + 1

    # Initialize score matrix and traceback matrix
    score_matrix = np.zeros((len_seq2, len_seq1))
    traceback_matrix = np.zeros((len_seq2, len_seq1, 2), dtype=int)

    max_score = 0
    max_pos = None

    # Fill in the matrices
    for i in range(1, len_seq2):
        for j in range(1, len_seq1):
            # Calculate the score for different options (match, mismatch, gap)
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[j - 1] == seq2[i - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            # Set the maximum score and update the traceback matrix
            score_matrix[i][j] = max(0, match, delete, insert)
            if score_matrix[i][j] == 0:
                traceback_matrix[i][j] = [0, 0]
            if score_matrix[i][j] == match:
                traceback_matrix[i][j] = [-1, -1]
            elif score_matrix[i][j] == delete:
                traceback_matrix[i][j] = [-1, 0]
            elif score_matrix[i][j] == insert:
                traceback_matrix[i][j] = [0, -1]

            # Track the maximum score and its position
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback to get the alignment
    if max_pos is not None:
        align1 = ""
        align2 = ""
        i, j = max_pos
        while i > 0 and j > 0 and score_matrix[i][j] > 0:
            if traceback_matrix[i][j][0] == -1 and traceback_matrix[i][j][1] == -1:
                align1 = seq1[j - 1] + align1
                align2 = seq2[i - 1] + align2
                i -= 1
                j -= 1
            elif traceback_matrix[i][j][0] == -1 and traceback_matrix[i][j][1] == 0:
                align1 = "-" + align1
                align2 = seq2[i - 1] + align2
                i -= 1
            elif traceback_matrix[i][j][0] == 0 and traceback_matrix[i][j][1] == -1:
                align1 = seq1[j - 1] + align1
                align2 = "-" + align2
                j -= 1

        return align1, align2, max_score
    else:
        return "", "", 0

def visualize_matrix(matrix):
    plt.figure(figsize=(8, 6))
    plt.imshow(matrix, cmap='Blues', interpolation='nearest')
    plt.colorbar()
    plt.xlabel('Sequence 1')
    plt.ylabel('Sequence 2')
    plt.title('Smith-Waterman Matrix')
    plt.show()

# Example usage:
seq1 = "AGTACGCA"
seq2 = "TATGC"
alignment1, alignment2, score = smith_waterman(seq1, seq2)
print("Sequence 1:", alignment1)
print("Sequence 2:", alignment2)
print("Score:", score)

# Visualization of the score matrix
score_matrix = np.zeros((len(seq2) + 1, len(seq1) + 1))
for i in range(len(seq2) + 1):
    for j in range(len(seq1) + 1):
        score_matrix[i][j] = smith_waterman(seq1[:j], seq2[:i])[2]

visualize_matrix(score_matrix)