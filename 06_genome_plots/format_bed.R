library(tidyverse)

"../03_alignment/rDNA_roi.bed" %>%
	read_tsv(col_names = F, col_types = "ciicic") %>%
	mutate(X6 = ".", X7 = X2, X8 = X3) %>%
	mutate(X9 = case_when(
		str_detect(X4, "(ETS|ITS)") ~ "105,226,68",
		str_detect(X4, "upstream_flank") ~ "234,51,35",
		TRUE ~ "44,102,61"
	)) %>%
	mutate(X4 = str_remove_all(X4, "::.+$")) %>%
	write_tsv("rDNA_roi.bed", col_names = F)
