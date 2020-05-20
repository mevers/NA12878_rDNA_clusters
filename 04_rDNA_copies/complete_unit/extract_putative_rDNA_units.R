library(tidyverse)


# Read FASTA index file to get lengths of contigs; when estimating the length
# of the last rDNA unit within a cluster based on the average length of the
# preceding rDNA units we need to make sure that we stay within the bounds of
# the contig
df <- "../../00_ref_sequences/assembly/albacore_canu_wtdbg_nanopolish2.fasta.fai" %>%
	read_tsv(col_names = FALSE, col_types = "cn---") %>%
	setNames(c("X1", "length"))


# Read BED files and find the loci where the first 100 bp rDNA fragment
# maps to; we will let these loci define the 5' start of a rDNA unit
bed <- list.files(
	path = "../../03_alignment",
	pattern = "^rDNA.+\\.bed",
	full.names = TRUE)
id <- bed %>% basename() %>% str_remove_all("\\.bed")
lst <- bed %>%
	setNames(id) %>%
	imap(~.x %>%
		read_tsv(col_names = FALSE, col_types = "ciicic") %>%
		filter(str_detect(X4, "(pos0|5'ETS)")) %>%
		left_join(df, by = "X1") %>%
		group_by(X1) %>%
		mutate(
			start = if_else(X6 == "+", X2, lag(X3)),
			end = if_else(X6 == "+", lead(X2), X3),
			name = "rDNA_start",
			value = X5) %>%
		mutate(
			end = ifelse(
				is.na(end),
				pmin(length, start + as.integer(mean(end - start, na.rm = TRUE))),
				end),
			start = ifelse(
				is.na(start),
				pmax(0, end - as.integer(mean(end - start, na.rm = TRUE))),
				start)) %>%
		select(X1, start, end, name, value, X6) %>%
	 	write_tsv(
			sprintf("rDNA_start_%s.bed", .y),
			col_names = FALSE,
			quote_escape = FALSE))
