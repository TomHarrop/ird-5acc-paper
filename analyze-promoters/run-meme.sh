#!/bin/bash


# last 100, both strands

# ./meme/bin/meme \
#        ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS.fasta \
#        -dna \
#        -revcomp \
#        -nmotifs 10\
#        -p 4 \
#        -o ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS-memeOUT



# last 100 forward strand

# ./meme/bin/meme \
#        ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS.fasta \
#        -dna \
#        -nmotifs 10\
#        -p 4 \
#        -o ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS-memeOUT-norev
       # -revcomp \


# top 1o0 forward strand

# ./meme/bin/meme \
#        ~/Desktop/ird-5acc-paper/seq/top_pc5_2000TSS.fasta \
#        -dna \
#        -nmotifs 10\
#        -p 4 \
#        -o ~/Desktop/ird-5acc-paper/seq/top_pc5_2000TSS-memeOUT-norev
       # -revcomp \


# last 100 forward strand with negative dataset

# ./meme/bin/meme \
#        ~/Desktop/ird-5acc-paper/seq/top_pc5_2000TSS.fasta \
#        -dna \
#        -nmotifs 10\
#        -p 4 \
#        -neg ~/Desktop/ird-5acc-paper/seq/random_genes_2000TSS.fasta \
#        -o ~/Desktop/ird-5acc-paper/seq/top_pc5_2000TSS-memeOUT-norev-neg
       # -revcomp \

# glam2 - gapped motives, not very powerful?

# ./meme/bin/glam2 \
#   -M n \
#   -o ~/Desktop/ird-5acc-paper/seq/top_pc5_2000TSS-glam2OUT \
#   ~/Desktop/ird-5acc-paper/seq/top_pc5_2000TSS.fasta

# top 400 forward strand

# ./meme/bin/meme \
#        ~/Desktop/ird-5acc-paper/seq/top_pc5_400_2000TSS.fasta \
#        -dna \
#        -nmotifs 20\
#        -p 4 \
#        -o ~/Desktop/ird-5acc-paper/seq/top_pc5_400_2000TSS-memeOUT-norev
# 
# ./meme/bin/meme \
#        ~/Desktop/ird-5acc-paper/seq/last_pc5_400_2000TSS.fasta \
#        -dna \
#        -nmotifs 20\
#        -p 4 \
#        -o ~/Desktop/ird-5acc-paper/seq/last_pc5_400_2000TSS-memeOUT-norev
#        
# send output to tomotom

# ./meme/bin/tomtom \
#     -oc ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS_tomtomOUT \
#     ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS-memeOUT-norev/meme.txt \
#     ~/Desktop/ird-5acc-paper/seq/JASPAR2018_CORE_plants_non-redundant.meme

# Mast

# ./meme/bin/mast \
#      -oc ~/Desktop/ird-5acc-paper/seq/last-pc5-mast-out \
#      ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS-memeOUT-norev/meme.txt \
#      ~/Desktop/ird-5acc-paper/seq/last_pc5_400_2000TSS.fasta
#      
# ./meme/bin/mast \
#      -oc ~/Desktop/ird-5acc-paper/seq/top-pc5-mast-out \
#      ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS-memeOUT-norev/meme.txt \
#      ~/Desktop/ird-5acc-paper/seq/top_pc5_400_2000TSS.fasta
#      
# ./meme/bin/mast \
#      -oc ~/Desktop/ird-5acc-paper/seq/random-mast-out \
#      ~/Desktop/ird-5acc-paper/seq/last_pc5_2000TSS-memeOUT-norev/meme.txt \
#      ~/Desktop/ird-5acc-paper/seq/random_genes_2000TSS.fasta