# seqalign

seqalign is a Rust library with tools for carrying out fast pairwise alignments of (for now) protein
sequences. At its core, seqalign implements the Smith-Waterman algorithm for sequence alignment,
with Gotoh's extension enabling affine gap penalties. This will eventually implement striped
Smith-Waterman (Farrar 2007), which greatly increases performance by employing SIMD vectorization to
calculate independent cells of the DP arrays in parallel.

This package comes with a CLI tool to enable pairwise alignment of one or more pairs of sequences as
part of whatever bioinformatics workflow the user is working on.

Plans for future development:
1. Implement global (i.e. Needleman-Wunsch) as well as local alignment.
2. Enable DNA alignment.
3. Parallelize across threads in addition to vector-level (SIMD).
