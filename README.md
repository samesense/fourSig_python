# fourSig_python
Incomplete port of fourSig to python.

I wanted to address the window size being fragment count, instead of
genomic size. Overlapping windows are used. Shuffling occurs by moving
fragment read counts, not reads.
