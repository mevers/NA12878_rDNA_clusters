#!/bin/bash

# conda install -c bioconda clustalo

for fn in ../../04_rDNA_copies/regions/*.fa; do
	bn=$(basename $fn)

	id=${bn/putative_rDNA_/}
	id=${id/.fa/}
	echo "Processing $id"

    clustalo \
    	-i $fn \
    	-t DNA \
    	-o msa_${bn%fa}clustal \
    	--outfmt clu \
    	--resno \
    	--wrap 60 \
    	-v

	params=(
		"-in clustal"
		"-conservation on"
		"-width 240"
		"-title \"Sequence similarity ${id} (relative to U13369.1 reference)\""
		"-bold"
		"-coloring mismatch"
		"-html head"
		"-css on"
		"-moltype dna"
		"msa_${bn%fa}clustal"
		"> msa_${bn%fa}html")

	# Highlight the top 10 uCNE elements in the LSU rRNA (28S)
	# Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/
	if [[ $id == "28S" ]]; then
		params+=(
			"-find \"CCGATAG:CCTAAG:CGTACC:TAACTT:GACTGTTTA:AAGACCC:TGGGGC:GGATAAC:GAGCTGGGTTTA:AGTACGAGAGGAAC\""
		)
	# Highlight the fwd/rev RT-qPCR primers for 5'ETS and ITS1
	# Reference: Rita Ferreira
	elif [[ $id == "5pETS" ]]; then
		params+=(
			"-find \"GCTCTTCGATCGATGTGGTGACG|GTCCTTCTCGCTCCGCCCG\""
		)
	elif [[ $id == "ITS1" ]]; then
		params+=(
			"-find \"GAGAACTCGGGAGGGAGAC|GAGAGAAAGAAGGGCGTGTC\""
		)
	fi

	eval ~/bin/mview "${params[@]}"

done
