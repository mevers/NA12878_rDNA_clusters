#!/bin/bash

# conda install -c bioconda pygenometracks

regions=(
	"tig00008411:1-310,000"
	"tig00008216:1-150,000"
	"tig00008215:1-154,000"
	"tig00002380:1-215,000")
ext=("pdf" "png")

for reg in ${regions[@]}; do
	fn=${reg/:/_}
	fn=${fn//,/}
	for e in ${ext[@]}; do
		pygenometracks \
			--tracks tracks.ini \
			--region $reg \
			--dpi 150 \
			--outFileName $fn.$e
	done
done
