#!/bin/bash

bowtie2-build \
	--threads 4 \
	../00_ref_sequences/assembly/albacore_canu_wtdbg_nanopolish2.fasta \
	albacore_canu_wtdbg_nanopolish2
