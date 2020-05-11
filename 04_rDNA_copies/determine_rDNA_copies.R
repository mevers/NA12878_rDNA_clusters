library(Rsamtools)
library(tidyverse)

bam <- list.files(
	path = "../03_alignment",
	pattern = "^rDNA.+\\.bam$",
	full.names = TRUE)
id <- basename(bam) %>% str_remove_all("(unsorted_|\\.bam)")

df <- setNames(bam, id) %>%
	map(~BamFile(.x) %>%
		scanBam() %>%
		pluck(1) %>%
		keep(names(.) %in% c("qname", "rname", "pos", "qwidth")) %>%
		as_tibble()) %>%
	bind_rows(.id = "Data")


# Plot number of hits per contig
df %>%
	group_by(Data, rname) %>%
	summarise(n_hits = n()) %>%
	ggplot(aes(rname, n_hits)) +
	geom_col() +
	facet_wrap(~ Data) +
	theme_minimal() +
	scale_y_continuous(
		sec.axis = sec_axis(
			~ . * 500 / 43000,
			name = "Rough rDNA copy number from crude inference")) + 
	labs(x = "Contig", y = "Number of hits") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("rDNA_frag_hits.pdf", width = 10, height = 6)



# Read BED files and find the loci where the first 100 bp rDNA fragment
# maps to; we will let these loci define the 5' start of a rDNA unit
bed <- list.files(
	path = "../03_alignment",
	pattern = "^rDNA.+\\.bed",
	full.names = TRUE)
id <- bed %>% basename() %>% str_remove_all("\\.bed")
lst <- bed %>%
	setNames(id) %>%
	imap(~.x %>%
		read_tsv(col_names = FALSE, col_types = "ciicic") %>%
		filter(str_detect(X4, "pos0")) %>%
		group_by(X1) %>%
		mutate(
			start = if_else(X6 == "+", X2, lag(X3)),
			end = if_else(X6 == "+", lead(X2), X3),
			name = "rDNA_start",
			value = X5) %>%
		select(X1, start, end, name, value, X6) %>%
		mutate(
			end = ifelse(
				is.na(end),
				start + as.integer(mean(end - start, na.rm = TRUE)),
				end),
			start = ifelse(
				is.na(start),
				pmax(0, end - as.integer(mean(end - start, na.rm = TRUE))),
				start)) %>%
	 	write_tsv(
			sprintf("rDNA_start_%s.bed", .y),
			col_names = FALSE,
			quote_escape = FALSE))
