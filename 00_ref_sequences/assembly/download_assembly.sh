#!/bin/bash

# Download FASTA file
curl -O https://gembox.cbcb.umd.edu/triobinning/albacore_canu_wtdbg_nanopolish2.fasta


# Index FASTA
samtools faidx albacore_canu_wtdbg_nanopolish2.fasta
