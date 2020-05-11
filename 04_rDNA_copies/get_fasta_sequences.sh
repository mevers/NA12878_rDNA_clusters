#!/bin/bash

bedtools getfasta \
	-s \
	-fi ../00_ref_sequences/assembly/albacore_canu_wtdbg_nanopolish2.fasta \
	-bed rDNA_start_rDNA_frags_len500_step500.bed \
	-fo putative_rDNA_seq_frags_len500_step500.fa

bedtools getfasta \
	-s \
	-fi ../00_ref_sequences/assembly/albacore_canu_wtdbg_nanopolish2.fasta \
	-bed rDNA_start_rDNA_frags_len1000_step1000.bed \
	-fo putative_rDNA_seq_frags_len1000_step1000.fa
