library("dplyr")
library("tidyr")
library("MatrixEQTL")
source("../MatrixEQTL/MatrixEQTL_wrapper.R")

output_file <- "---"
phenotype_file <- "---"
#q16_pca.R abil toodetud fail
covariates_file <- "---"
#q16_selected.R abil toodetud fail
sample_file <- "---"
#vcfist_snpspos.txt abil toodetud fail
snpspos_file <- "---"
#vcfist_doosimaatriks.txt abil toodetud fail
genotypes_file <- "---"
#q16_expression_chunk.R abil toodetud fail
expression_file <- "---"

kept_samples <- read.table(sample_file, sep="\t", header=TRUE)
kept_samples_by_genotype_id <- kept_samples[order(kept_samples$genotype_id),]

phenotypes <- read.table(phenotype_file, sep="\t", header = TRUE)

snpspos <- read.table(snpspos_file, sep="\t", header=TRUE, row.names = NULL)
colnames(snpspos) <- c("snpid", "chr", "pos")
snpspos$chr <- as.character(snpspos$chr)

genotypes = read.table(genotypes_file, sep="\t", header=TRUE)
geno_names = paste("chr",genotypes$CHROM, "_", genotypes$POS, "_", genotypes$REF, "_", genotypes$ALT, sep="")
genotypes = genotypes[,-(1:4)]
substr(colnames(genotypes), 4, 4) <- "@"

geno_data <- genotypes[,colnames(genotypes) %in% kept_samples$genotype_id]
#sorteerin tulbad tähestiku järjekorras
geno_data <- geno_data %>% 
  dplyr::select(sort(names(.)))
geno_data <- as.matrix(geno_data)
colnames(geno_data) <- kept_samples_by_genotype_id$genotype_id
rownames(geno_data) <- geno_names

expression <- read.table(expression_file, sep="\t", header=TRUE, row.names=NULL)

#jätan phenotype_id välja
exp_data <- expression[,2:length(expression)]
#sorteerin tulbad tähestiku järjekorras
exp_data <- exp_data %>% 
  dplyr::select(sort(names(.)))

#Yang et al, inverse normal transformation
l <- length(colnames(exp_data))
m <- length(rownames(exp_data))
df <- data.frame(matrix(NA, nrow=m, ncol=1))

exp_data <- t(as.matrix(exp_data))
exp_data <- qnorm((rank(exp_data, na.last = "keep") - 0.5) / sum(!is.na(exp_data)))


for (colnr in 1:l) {
  indices <- seq(colnr,length(exp_data),l)
  colname <- kept_samples_by_genotype_id$genotype_id[colnr]
  df[[colname]] <- exp_data[indices]
}

df <- df[,-1]

exp_data <- as.matrix(df)

rownames(exp_data) <- expression$phenotype_id
colnames(exp_data) <- kept_samples_by_genotype_id$genotype_id


gene_pos <- data_frame(geneid=phenotypes$gene_id, chr=phenotypes$chromosome, left=phenotypes$gene_start, right=phenotypes$gene_end)

#pca.R abil toodetud fail
covariates <- read.table(covariates_file, sep="\t", header=TRUE)
covariates_cols <- colnames(covariates)
covariates <- as.matrix(covariates)
colnames(covariates) <- covariates_cols
substr(colnames(covariates), 4, 4) <- "@"

runMatrixEQTL(exp_data, geno_data, snpspos, gene_pos, covariates=covariates, cisDist = 5e5, outputFileName = output_file, pvOutputThreshold = 1, permute = FALSE, model = modelLINEAR)
