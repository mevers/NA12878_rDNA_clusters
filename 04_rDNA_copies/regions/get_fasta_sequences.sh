#!/bin/bash

assembly="../../00_ref_sequences/assembly/albacore_canu_wtdbg_nanopolish2.fasta"
rDNA="../../02_rDNA_frags/rDNA_transcribed.fa.gz"

for fn in *.bed; do
	id=${fn/putative_rDNA_/}
	id=${id/.bed/}
	id=${id/5p/5\'}
	id=${id/3p/3\'}

	zcat < $rDNA | grep $id -A 1 > tmp1.fa
	bedtools getfasta \
		-s \
		-fi $assembly \
		-bed $fn \
		-fo tmp2.fa
	cat tmp1.fa tmp2.fa > ${fn%bed}fa
done

# Clean up
rm -f tmp1.fa
rm -f tmp2.fa
