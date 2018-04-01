# RecNW: Recursive Needleman-Wunsch aligner

RecNW is a Python package written in C++ and compiled with Cython that can be used either as a Python command line tool or a Python function. this repository contains the code for the paper, "RecNW: A fast recurrent pairwise aligner for targeted sequencing", by Alexandre Yahi, Tuuli Lappalainen, Pejman Mohammadi*, and Nicholas P Tatonetti* (co-senior authors).
RecNW is particularly useful for analyzing targeted sequencing data where reads need to be aligned to a specific reference sequence. RecNW leverages reads similarities by: 1/ deduplicating the reads and aligning only unique reads, 2/ sorting reads and recursively re-using blocs of the alignment matrices to speed up computation.

## Quick start

Primary dependencies: 


### Compiling the package with Cython


### Aligning the sample datasets

## Content of the package

The RecNW `recnw.so` Python package contains four functions: `nw_aff`, `nw_lin`, `recnw_aff`, `recnw_lin`. All functions offer the possibility to make head and tail gaps penalty-free for both reads independently.

### `nw_aff`

This function is an implementation of the Needleman-Wunsch algorithm with affine gap penalty (Gotoh algorithm).

## About the algorithms

Let assume we have two sequences, seq1 and seq2 of respective lengths n and m.

### Linear gap penalty


### Affine gap penalty
