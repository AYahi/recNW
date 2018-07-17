# RecNW: Recursive Needleman-Wunsch aligner

RecNW is a Python package written in C++ and compiled with Cython that can be used either as a Python command line tool or a Python function. this repository contains the code for the paper, _[RecNW: A fast recurrent pairwise aligner for targeted sequencing]()_, by Alexandre Yahi, Tuuli Lappalainen, Pejman Mohammadi<sup>+</sup>, and Nicholas P Tatonetti<sup>+</sup> (<sup>+</sup>co-senior authors).

RecNW is particularly useful for analyzing targeted sequencing data where reads need to be aligned to a specific reference sequence. RecNW leverages reads similarities by: 1/ deduplicating the reads and aligning only unique reads, 2/ sorting reads and recursively re-using blocs of the alignment matrices to speed up computation.

## Quick start

Primary dependencies: `Cython`.
The sample code and command line tool were written for Python 2.7.

To install RecNW and compile it for your machine (Linux/Mac):

```
git clone git@github.com:AYahi/recNW.git
cd recNW
cd cython_recnw
python setup.py build_ext --inplace
```

### Aligning the sample datasets

To run RecNW on a sample dataset provided in the repository with affine gap penalty and default parameters with the command line version:
```
python recnw.py --input_file data/sim_sample_1.fastq.gz --ref data/ACCS.fasta --penalty 'aff'
```
The results can be found in the `/results` directory and have the following formatting:
```
>seq1_1
Score:73.0
GCTCTGCCAGGCCCTAGTGCGC
||||||||||||||||||||||
GCTCT-----GCCCTAGTGCGC
>seq1_2
Score:52.0
GCTCTGCCAGGCCCTAGTGCGC
||||||||||||||||||||||
GCTCTGC---GACCTTGTG---
...
```

Please note that the Python function when used inside a Python script returns the results in a list of tuples respecting the order of the input sequences.

#### Parameters

|Input parameters|Default|Description|
|----------------|----|-----------|
| --input_file | None | Path to the `.fasta` or `.fastq` file with the reads to align. `gzip` supported.|
| --output_dir | 'results' | Directory for the output.|
| --ref | None | Path to the `.fasta` reference file, or full reference sequence.|
| --penalty | `aff` | `aff` for affine gap penalty, `lin` for linear gap penalty.|
| --gap_op | 8 | Gap opening penalty, if affine gap penalty. Must be positive.|
| --gap_ext | 1 | Gap extension penalty, if affine gap penalty. Must be positive.|
| --gap_penalty | 8 | Gap penalty, if linear gap penalty. Must be positive.|
| --match | 5 | Score of a match between nucleotides.|
| --mismatch | -4 | Score of a mismatch between nucleotides.|
| --head_free | (False, True)| Free head gap penalty for reference, and read to align, respectively.|
| --tail_gap | (False, True)| Free tail gap penalty for reference, and read to align, respectively.|
| --memory | True | recNW version re-using blocs of the alignment matrices for speed increase with optimal alignment. False for classic Needleman-Wunsch.|



## Content of the Python package

The RecNW `recnw.so` Python package contains four functions: `nw_aff`, `nw_lin`, `recnw_aff`, `recnw_lin`. All functions offer the possibility to make head and tail gaps penalty-free for both reads independently.

### Aligning a single read
### `nw_aff`

This function is an implementation of the Needleman-Wunsch algorithm with affine gap penalty (Gotoh algorithm). It aligns **a single read** to a given reference and return their alignments and score.

|Input parameters|Type|Description|
|----------------|----|-----------|
|ref| string | Reference sequence against the read must be aligned|
|seq| string | Sequence from the reads dataset to align with the reference|
|gap_op (default:8)| float | Gap opening penalty, must be > 0. (affine gap penalty)|
|gap_ext (default:1)| float | Gap extension penalty, must be > 0 (affine gap penalty)|
|match (default:5)| float | Score of a match between two symbols, usually positive|
|mismatch (default:-4)| float | Score of a mismatch, usually negative|
|head_free (default: (False, True))|(bool, bool)| Option to have free head gap penalty for the reference, and the sequence to align.|
|tail_free (default: (False, True))|(bool, bool)| Option to have free tail gap penalty for the reference, and the sequence to align.|
|sim (default:-1)| int | -1 to initialize the alignment matrices, indicates the number of similar nucleotides with the previous aligned sequence (to use when memory is kept)|
|terminate (default:1)| int | 1 to delete matrices and de-allocate memory, 0 to keep the alignment matrices when used with the `sim` parameter.|

|Output|Type|Description|
|------|----|-----------|
|al1|string|alignment of the reference, with '-' for gaps|
|al2|string|alignment of the read sequence, with '-' for gaps|
|score|float| score of the alignment|

#### `nw_lin`

This function is an implementation of the Needleman-Wunsch algorithm with linear gap penalty. It aligns **a single read** to a given reference and return their alignments and score.

|Input parameters|Type|Description|
|----------------|----|-----------|
|ref| string | Reference sequence against the read must be aligned|
|seq| string | Sequence from the reads dataset to align with the reference|
|gap_penalty (default:8)| float | Gap penalty, must be > 0.|
|match (default:5)| float | Score of a match between two symbols, usually positive|
|mismatch (default:-4)| float | Score of a mismatch, usually negative|
|head_free (default: (False, True))|(bool, bool)| Option to have free head gap penalty for the reference, and the sequence to align.|
|tail_free (default: (False, True))|(bool, bool)| Option to have free tail gap penalty for the reference, and the sequence to align.|
|sim (default:-1)| int | -1 to initialize the alignment matrices, indicates the number of similar nucleotides with the previous aligned sequence (to use when memory is kept)|
|terminate (default:1)| int | 1 to delete matrices and de-allocate memory, 0 to keep the alignment matrices when used with the `sim` parameter.|

|Output|Type|Description|
|------|----|-----------|
|al1|string|alignment of the reference, with '-' for gaps|
|al2|string|alignment of the read sequence, with '-' for gaps|
|score|float| score of the alignment|

### Aligning a dataset

Supported dataset formats include `.fasta` and `.fastq` files. The files can be compressed with `gzip`.
These functions respectively implement an affine gap penalty strategy, and a linear gap penalty. They deduplicate the reads for alignment purposes, and the user has the option to use the `recNW` option to keep alignment matrices in memory and re-use relevant blocks of these matrices based on the reads similarity with the `memory` parameter.

#### `recnw_aff`

This function uses `nw_aff` to perform the Needleman-Wunsch alignment of **a whole dataset** against a given reference with affine gap penalty (Gotoh algorithm).

|Input parameters|Type|Description|
|----------------|----|-----------|
|file_path| string | Relative path to the file to align |
|seq| string | Relative path to the fasta file with the reference, or reference sequence string|
|output_dir| string | Relative path to the directory where the alignment text file will be saved.|
|gap_op (default:8)| float | Gap opening penalty, must be > 0. (affine gap penalty)|
|gap_ext (default:1)| float | Gap extension penalty, must be > 0 (affine gap penalty)|
|match (default:5)| float | Score of a match between two symbols, usually positive|
|mismatch (default:-4)| float | Score of a mismatch, usually negative|
|head_free (default: (False, True))|(bool, bool)| Option to have free head gap penalty for the reference, and the sequence to align.|
|tail_free (default: (False, True))|(bool, bool)| Option to have free tail gap penalty for the reference, and the sequence to align.|
|memory (default: True)| bool | Use recNW recursive alignment matrix bloc re-use|

**Output** List of tuples (alignment of the sequence, alignment of the reference, score, sequence name)

#### `recnw_lin`

This function uses `nw_lin` to perform the Needleman-Wunsch alignment of **a whole dataset** against a given reference with linear gap penalty.

|Input parameters|Type|Description|
|----------------|----|-----------|
|file_path| string | Relative path to the file to align |
|seq| string | Relative path to the fasta file with the reference, or reference sequence string|
|output_dir| string | Relative path to the directory where the alignment text file will be saved.|
|gap_penalty (default:8)| float | Gap penalty, must be > 0.|
|match (default:5)| float | Score of a match between two symbols, usually positive|
|mismatch (default:-4)| float | Score of a mismatch, usually negative|
|head_free (default: (False, True))|(bool, bool)| Option to have free head gap penalty for the reference, and the sequence to align.|
|tail_free (default: (False, True))|(bool, bool)| Option to have free tail gap penalty for the reference, and the sequence to align.|
|memory (default: True)| bool | Use recNW recursive alignment matrix bloc re-use|

**Output** List of tuples (alignment of the sequence, alignment of the reference, score, sequence name)
