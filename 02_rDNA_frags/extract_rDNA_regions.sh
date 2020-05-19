#!/bin/bash

id="U13369.1"
flank=500


# Extract coding regions from annotation
cat ../00_ref_sequences/rDNA_GenBank/${id}.bed | grep \
    -e "\dS" -e "ETS" -e "ITS" > ${id}_transcribed.bed


# Create upstream flank region from last line of BED file and store in
# final BED file for the regions of interest (ROI)
function join_by { local IFS="$1"; shift; echo "$*"; }
last=$(tail -n 1 ../00_ref_sequences/rDNA_GenBank/${id}.bed)
arr=(${last//\t/})
arr[1]=$(( ${arr[2]} - $flank ))
arr[3]="upstream_flank"
arr[6]=${arr[1]}
join_by $'\t' "${arr[@]}" > tmp.bed
cat ${id}_transcribed.bed tmp.bed > ${id}_roi.bed
rm -f tmp.bed


# Get sequences for ROI's
bedtools getfasta \
    -name \
	-fi ../00_ref_sequences/rDNA_GenBank/${id}.fa \
	-bed ${id}_roi.bed | bgzip > rDNA_roi.fa.gz
