from align_sequence import globa_alignment, affine_gap_alignment, local_alignment
from align_sequence import print_global_alignment, print_affine_gap_alignment, print_local_alignment,print_alignment_matrix
from read_fasta import read_file
import time

import os
import psutil

pid = os.getpid()
py = psutil.Process(pid)

sequences = read_file("sequences.txt")

affine_gap_g1_seq1 = sequences[0]
affine_gap_g1_seq2 = sequences[1]
affine_gap_g2_seq1 = sequences[2]
affine_gap_g2_seq2 = sequences[3]
local_g1_seq1 = sequences[4]
local_g1_seq2 = sequences[5]
local_g2_seq1 = sequences[6]
local_g2_seq2 = sequences[7]

def compare_affine_gap_alignment_to_global_alignment(seq1,seq2):
    global_align = globa_alignment(seq1, seq2)
    top,mid,bot = affine_gap_alignment(seq1,seq2)
    print_global_alignment(global_align,seq1,seq2)
    print("Global alignment score of ", global_align[-1][-1])
    print_affine_gap_alignment(top,mid,bot,seq1,seq2)
    print("Affine gap alignment score of ", mid[-1][-1])

def compare_local_alignment_to_global_alignment(seq1,seq2):
    global_align = globa_alignment(seq1, seq2)
    local = local_alignment(seq1,seq2)
    print_global_alignment(global_align,seq1,seq2)
    print("Global alignment score of ", global_align[-1][-1])
    print_local_alignment(local,seq1,seq2)
    print("Local alignment score of ", local[-1][-1])

compare_affine_gap_alignment_to_global_alignment(affine_gap_g1_seq1,affine_gap_g1_seq2)
compare_affine_gap_alignment_to_global_alignment(affine_gap_g2_seq1,affine_gap_g2_seq2)
compare_local_alignment_to_global_alignment(local_g1_seq1,local_g1_seq2)
compare_local_alignment_to_global_alignment(local_g2_seq1,local_g2_seq2)