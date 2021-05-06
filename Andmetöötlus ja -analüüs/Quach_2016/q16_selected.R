#q16_add_PC1_to_metadata.R abil toodetud fail
meta_file_with_PC1 <- "---"
output_file_pop1 <- "---"
output_file_pop2 <- "---"
quach <- read.table(meta_file_with_PC1, sep="\t", header=TRUE)

LPS_bool <- quach$condition == "LPS"
rna_bool <- quach$rna_qc_passed
genotype_bool <- quach$genotype_qc_passed
population_bool <- quach$geno_PC1 < 0


population_1 <- quach[LPS_bool & rna_bool & genotype_bool & population_bool,]
population_2 <- quach[LPS_bool & rna_bool & genotype_bool & !population_bool,]

write.table(output_file_pop1, file="population_1.tsv", sep="\t", quote = FALSE, row.names = FALSE)
write.table(output_file_pop2, file="population_2.tsv", sep="\t", quote = FALSE, row.names = FALSE)