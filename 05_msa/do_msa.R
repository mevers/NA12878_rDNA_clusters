library(msa)

fn <- "../04_rDNA_copies/putative_rDNA_seq_frags_len500_step500.fa"
fa <- readDNAStringSet(fn)

aln <- msa(fa)
save(aln, file = "msa_putative_rDNA_seq_frags_len500_step500.RData")

chunk_size <- 2090
for (start in seq(1, ncol(aln), by = chunk_size)) {
	end <- min(start + chunk_size - 1, ncol(aln))
	alnPart <- DNAMultipleAlignment(subseq(unmasked(aln), start, end))
	consensus <- msaConsensusSequence(alnPart)
	msaPrettyPrint(
		alnPart,
		output = "pdf",
		subset = NULL,
		showLogo = "top",
		showNumbering = "none",
		askForOverwrite = FALSE,
		showLegend = FALSE,
		consensusThreshold = c(50, 80),
		consensusColors = "Gray",
		furtherCode=c(
			"\\defconsensus{.}{lower}{upper}",
			sprintf("\\showruler{top}{%i}", start)),
		file = sprintf("aln_pos%05i-%05i.pdf", start, end))
}



msaPrettyPrint(
	alnPart,
	output = "pdf",
	subset = NULL,
	showLogo = "top",
	#showNumbering = "none",
	askForOverwrite = FALSE,
	showLegend = FALSE,
	consensusThreshold = c(50, 80),
	consensusColors = "Gray",
	furtherCode=c(
		"\\defconsensus{.}{lower}{upper}",
		"\\showruler{100}{top}"),
	file = "tmp.pdf")
