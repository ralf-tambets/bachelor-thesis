library("dplyr")

#meta failile PC1 lisamiseks ja populatsioonide eristamiseks

PCA_file <- "---"
output_file <- "---"
meta_file <- "---"
PCA <- read.table(PCA_file, sep="\t", header=TRUE)
n <- PCA$genotype_id
m <- colnames(PCA)[-1]
substr(m, 4, 4) <- "@"

PCA <- as.data.frame(t(PCA[,-1]))
colnames(PCA) <- n

joiner <- data.frame(genotype_id=m,geno_PC1=PCA$geno_PC1)

metadata <- read.table(meta_file, sep="\t", header=TRUE)
metadata <- left_join(metadata, joiner, by="genotype_id", copy=FALSE)

write.table(metadata, output_file, sep="\t", quote=FALSE, row.names=FALSE)