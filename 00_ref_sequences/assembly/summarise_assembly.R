library(Biostrings)
library(tidyverse)

fa <- readDNAStringSet("albacore_canu_wtdbg_nanopolish2.fasta")

data.frame(size = width(fa) / 1000) %>%
	ggplot(aes(size)) +
	geom_histogram(bins = 100) +
	theme_minimal() +
	scale_x_log10() +
	labs(
		title = "NA12878-based Canu 1.7 + WTDBG + Nanopolish reference assembly",
		subtitle = sprintf("Number of contigs: %i", length(fa)),
		x = "Size of contig in kb", y = "Count")
ggsave("hist_contig_size.pdf", height = 6, width = 8)
ggsave("hist_contig_size.png", height = 6, width = 8)
