#!/bin/bash

# conda install -c bioconda clustalo
# Run (this takes around 40 mins):
clustalo \
	-i ../04_rDNA_copies/putative_rDNA_seq_frags_len500_step500.fa \
	-t DNA \
	-o msa_putative_rDNA.clustal \
	--outfmt clu \
	--resno \
	--wrap 60 \
	--threads 4 \
	-v


# https://desmid.github.io/mview/
# Run (this takes around 15 mins):
~/bin/mview \
	-in clustal \
	-conservation on \
	-width 120 \
	-title "Inferred rDNA unit sequence similarity" \
	-bold \
	-coloring group \
	-threshold 90 \
	-consensus on \
	-ignore class \
	-html head \
	-css on \
	-moltype dna \
	msa_putative_rDNA.clustal > msa_putative_rDNA.html


# MSA relative to U13369.1
cat \
	../00_ref_sequences/rDNA_GenBank/U13369.1.fa \
	../04_rDNA_copies/putative_rDNA_seq_frags_len500_step500.fa > all_rDNA.fa
clustalo \
	-i all_rDNA.fa \
	-t DNA \
	-o msa_all_rDNA.clustal \
	--outfmt clu \
	--resno \
	--wrap 60 \
	--threads 4 \
	-v


~/bin/mview \
	-in clustal \
	-conservation on \
	-width 120 \
	-title "Inferred rDNA unit sequence similarity (including U13369.1 reference)" \
	-bold \
	-coloring mismatch \
	-html head \
	-css on \
	-moltype dna \
	msa_all_rDNA.clustal > msa_all_rDNA.html
