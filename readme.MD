# Phylogenetic Tree using Smith Waterman Algorithm

## Smith Waterman Algorithm

The [SW.py](https://github.com/fl-sll/SWPhylogenetics/blob/main/SW.py) file contains hard coded python code containing the Smith Waterman algorithm which uses a dynamic programming approach. It calculates distances from 1 local alignment to another with match score, mismatch penalty, and gap penalty. More detailed explanation can be seen [here](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm). The code outputs a matrix showing the Smith Waterman scores.

## Phylogenetic Tree using BioPython aligner

The SW algorithm used in [PT.py](https://github.com/fl-sll/SWPhylogenetics/blob/main/PT.py) is made using the help of pairwise2.align.localxx from the BioPython module. Same match score, mismatch penalty, and gap penalty is being used as in [SW.py](https://github.com/fl-sll/SWPhylogenetics/blob/main/SW.py). We then plot the dendogram using heirarchical clustering using the UPGMA (Unweightied Pair-Group Method with Arithmetic mean) algorithm and the help of the Matplotlib module, we are able to plot the phylogenetic tree.