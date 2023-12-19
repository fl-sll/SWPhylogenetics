import numpy as np
import pandas as pd

def align_global(seqi, seqii, match, mismatch, open_gap, gap):
        
    ni = len(seqi)
    nii = len(seqii)

    # initialization of NW matrix

    M = np.zeros((ni + 1, nii + 1))
    M[0, 1] += open_gap
    M[1, 0] += open_gap

    M[1:, 0] = np.linspace(-2, -2 + (ni - 1) * gap, ni)
    M[0, 1:] = np.linspace(-2, -2 + (nii - 1) * gap, nii)

    # Compute the NeeW matrix

    for i in range(1, ni + 1):
        for j in range(1, nii + 1):
            diag = M[i - 1, j - 1]
            ver = M[i - 1, j]
            hor = M[i, j - 1]
            diag += match if (seqi[i - 1] == seqii[j - 1]) else mismatch
            ver += open_gap if (i == 1) else gap
            hor += open_gap if (j == 1) else gap
            M[i, j] = max([diag, ver, hor])
    

    # Find the optimal alignment

    al_seqi = []
    al_seqii = []
    al_seqi.append(seqi[-1])
    al_seqii.append(seqii[-1])


    i = ni - 1
    j = nii - 1


    while i > 0 and j > 0:
        diag = M[i - 1, j - 1]
        ver = M[i - 1, j]
        hor = M[i, j - 1]
        

        if diag >= ver and diag >= hor:
            i -= 1
            j -= 1
            al_seqi.append(seqi[i])
            al_seqii.append(seqii[j])

        elif hor > diag and hor > ver:
            j -= 1
            al_seqii.append(seqii[j])
            al_seqi.append('-')


        elif ver > diag and ver > hor:
            i -= 1
            al_seqi.append(seqi[i])
            al_seqii.append('-')
            
            
        
    al_seqi = ''.join(al_seqi)[::-1]
    al_seqii = ''.join(al_seqii)[::-1]
            

    return [al_seqi, al_seqii]


def align_local(seqi, seqii, match, mismatch, open_gap, gap):

    ni = len(seqi)
    nii = len(seqii)

    # Initialize the Smith-Waterman matrix

    M = np.zeros((ni + 1, nii + 1))

    # Compute the Smith-Waterman matrix

    for i in range(1, ni + 1):
        for j in range(1, nii + 1):
            diag = M[i - 1, j - 1]
            ver = M[i - 1, j]
            hor = M[i, j - 1]
            diag += match if (seqi[i - 1] == seqii[j - 1]) else mismatch
            ver += open_gap if (i == 1) else gap
            hor += open_gap if (j == 1) else gap
            M[i, j] = max([diag, ver, hor, 0])

    # Find the optimal local alignment
    i, j = np.unravel_index(M.argmax(), M.shape)
    al_seqi = []
    al_seqii = []

    while i > 0 and j > 0 and M[i, j] > 0:
        diag = M[i - 1, j - 1]
        ver = M[i - 1, j]
        hor = M[i, j - 1]
        

        if diag >= ver and diag >= hor:
            i -= 1
            j -= 1
            al_seqi.append(seqi[i])
            al_seqii.append(seqii[j])

        elif hor > diag and hor > ver:
            j -= 1
            al_seqii.append(seqii[j])
            al_seqi.append('-')


        elif ver > diag and ver > hor:
            i -= 1
            al_seqi.append(seqi[i])
            al_seqii.append('-')
            
            
    al_seqi = ''.join(al_seqi)[::-1]
    al_seqii = ''.join(al_seqii)[::-1]

    return [al_seqi, al_seqii]


def align(seqi, seqii, alignment_type = 'global', match = 2, mismatch = -1, open_gap = -2, gap = -1):

    seqi = np.array(list(seqi))
    seqii = np.array(list(seqii))


    # substitution matrix

    # matrix = [[2, -6, -6, -6], [-6, 2, -6, -6], 
    #           [-6, -6, 2, -6], [-6, -6, -6, 2]]

    if alignment_type == 'local':
        return align_local(seqi, seqii, match, mismatch, open_gap, gap)
    elif alignment_type == 'global':
        return align_global(seqi, seqii, match, mismatch, open_gap, gap)
    else:
        raise ValueError('Invalid alignment type: {}'.format(alignment_type))

# seq1 = 'CACACAGTGACTAGCTAGCTACGATC'
# seq2 = 'CACACAGTCGACTAGCTAGCACGATC'

sequences = [
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAACGAAAAAAAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAATTAAAAAAAAAAAAAAAAAAAAAAA"
]

num_seqs = len(sequences)
alignment = []
al_seqs = len(alignment)
aligns = []
aligned = []

# for i in range(num_seqs):
#     for j in range(i + 1, num_seqs):
#         alignment.append(align(sequences[i], sequences[j], 'global'))

for i in range(num_seqs):
    for j in range(i + 1, num_seqs):
        alignment.append(align(sequences[i], sequences[j], 'local'))

# print(alignment[0][0])

for i in alignment:
    aligns.append(i[1])
        
# ml = 0

# for i in aligns:
#     if len(i) > ml:
#         ml = len(i)

# for i in aligns:
#     if len(i) < ml:
#         aligned.append(i + "-" * (ml - len(i)))
    
ml = 0
for i in aligns:
    if len(i) > ml:
        ml = len(i)
    
for i in aligns:
    if len(i) < ml:
        aligned.append(i + "-" * (ml - len(i)))
    else:
        aligned.append(i)

aligned = list(dict.fromkeys(aligned))
    
for i in aligned:
    print(i)



# print(*align(seq1, seq2, 'global'), sep='\n')