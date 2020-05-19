library(tidyverse)

"../../03_alignment/rDNA_roi.bed" %>%
	read_tsv(col_names = FALSE, col_types = "ciicic") %>%
	mutate(
		X4 = str_remove_all(X4, "::.*$"),
		group = str_replace(X4, "'", "p")) %>%
	group_by(group) %>%
	nest() %>%
	mutate(tmp = map2(
		group, data,
		~write_tsv(.y, sprintf("putative_rDNA_%s.bed", .x), col_names = FALSE)))
