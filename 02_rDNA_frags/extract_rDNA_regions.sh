#!/bin/bash

id="U13369.1"


# Extract coding regions from annotation
cat ../00_ref_sequences/rDNA_GenBank/${id}.bed | grep \
    -e "\dS" -e "ETS" -e "ITS" > ${id}_transcribed.bed


# Get coding region sequences
bedtools getfasta \
    -name \
	-fi ../00_ref_sequences/rDNA_GenBank/${id}.fa \
	-bed ${id}_transcribed.bed | bgzip > rDNA_transcribed.fa.gz
