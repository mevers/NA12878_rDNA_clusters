library(Rsamtools)
library(tidyverse)
library(gridExtra)

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
df_sum <- df %>%
	group_by(Data, rname) %>%
	summarise(n_hits = n())
gg1 <- df_sum %>%
	filter(Data == id[1]) %>%
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
gg2 <- df_sum %>%
	filter(Data == id[2]) %>%
	ggplot(aes(rname, n_hits)) +
	geom_col() +
	facet_wrap(~ Data) +
	theme_minimal() +
	scale_y_continuous(
		sec.axis = sec_axis(
			~ . / 8,
			name = "Rough rDNA copy number from crude inference")) +
	labs(x = "Contig", y = "Number of hits") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gg <- grid.arrange(gg1, gg2, ncol = 2)
ggsave("rDNA_frag_hits.pdf", gg, width = 10, height = 6)
ggsave("rDNA_frag_hits.png", gg, width = 10, height = 6)
