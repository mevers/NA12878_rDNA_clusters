#!/bin/bash

bowtie2 \
	-x ../01_bowtie2_ref/albacore_canu_wtdbg_nanopolish2 \
	-U ../02_rDNA_frags/rDNA_frags_len500_step500.fa.gz \
	-f \
	--threads 4 \
	--all \
	2> ../logs/bowtie2_rDNA_frags_len500_step500.log | samtools view \
	-bSh -@ 4 -F4 - > unsorted_rDNA_frags_len500_step500.bam
bowtie2 \
	-x ../01_bowtie2_ref/albacore_canu_wtdbg_nanopolish2 \
	-U ../02_rDNA_frags/rDNA_transcribed.fa.gz \
	-f \
	--threads 4 \
	--all \
	2> ../logs/bowtie2_rDNA_transcribed.log | samtools view \
	-bSh -@ 4 -F4 - > unsorted_rDNA_transcribed.bam


samtools sort -@ 4 \
	-o rDNA_frags_len500_step500.bam unsorted_rDNA_frags_len500_step500.bam
samtools sort -@ 4 \
	-o rDNA_transcribed.bam unsorted_rDNA_transcribed.bam


samtools index rDNA_frags_len500_step500.bam
samtools index rDNA_transcribed.bam


bamCoverage \
	-b rDNA_frags_len500_step500.bam \
	-o rDNA_frags_len500_step500.bw
#bamCoverage \
#	--samFlagExclude 16 \
#	-b rDNA_frags_len500_step500.bam \
#	-o rDNA_frags_len500_step500_pos.bw
#bamCoverage \
#	--samFlagInclude 16 \
#	-b rDNA_frags_len500_step500.bam \
#	-o rDNA_frags_len500_step500_neg.bw
bamCoverage -b rDNA_transcribed.bam -o rDNA_transcribed.bw


bedtools bamtobed -i rDNA_frags_len500_step500.bam > rDNA_frags_len500_step500.bed
bedtools bamtobed -i rDNA_transcribed.bam > rDNA_transcribed.bed
