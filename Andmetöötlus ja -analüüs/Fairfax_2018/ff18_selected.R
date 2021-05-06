meta_file <- ""

fairfax <- read.table(meta_file, sep="\t", header=TRUE)

LPS24_bool <- fairfax$condition == "LPS24"
rna_bool <- fairfax$rna_qc_passed
genotype_bool <- fairfax$genotype_qc_passed

fairfax_selected <- fairfax[LPS24_bool & rna_bool & genotype_bool,]

write.table(fairfax_selected, file="ff18_selected.tsv", sep="\t", quote = FALSE, row.names = FALSE)
