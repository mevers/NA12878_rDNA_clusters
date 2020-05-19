#!/bin/bash

id=("rDNA_frags_len500_step500" "rDNA_roi")

for name in "${id[@]}"; do

	bowtie2 \
		-x ../01_bowtie2_ref/albacore_canu_wtdbg_nanopolish2 \
		-U ../02_rDNA_frags/$name.fa.gz \
		-f \
		--threads 4 \
		--all \
		2> ../logs/bowtie2_$name.log | samtools view \
		-bSh -@ 4 -F4 - > unsorted_$name.bam

	samtools sort -@ 4 \
		-o $name.bam unsorted_$name.bam

	rm -f unsorted_$name.bam

	samtools index $name.bam

	bedtools bamtobed -i $name.bam > $name.bed

done
