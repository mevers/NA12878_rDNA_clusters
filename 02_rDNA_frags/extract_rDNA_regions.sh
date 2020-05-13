#!/bin/bash

bedtools getfasta \
    -name \
	-fi ../00_ref_sequences/rDNA_GenBank/U13369.1.fa \
	-bed U13369.1_transcribed.bed | bgzip > rDNA_transcribed.fa.gz
