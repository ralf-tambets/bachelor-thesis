library("dplyr")
library("tidyr")
library("MatrixEQTL")
source("../MatrixEQTL/MatrixEQTL_wrapper.R")

output_file <- "---"
juhtvariant_file <- "---"
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

snpspos <- read.table(snpspos_file, sep="\t", header=TRUE, row.names = NULL)
colnames(snpspos) <- c("snpid", "chr", "pos")
snpspos$chr <- as.character(snpspos$chr)

genotypes = read.table(genotypes_file, sep="\t", header=TRUE)
geno_names = paste("chr",genotypes$CHROM, "_", genotypes$POS, "_", genotypes$REF, "_", genotypes$ALT, sep="")
substr(colnames(genotypes), 4, 4) <- "@"

geno_data <- genotypes[,colnames(genotypes) %in% kept_samples$genotype_id]
#sorteerin tulbad tähestiku järjekorras
geno_data <- geno_data %>% 
  dplyr::select(sort(names(.)))
geno_data <- as.matrix(geno_data)
colnames(geno_data) <- kept_samples_by_genotype_id$genotype_id
rownames(geno_data) <- geno_names

expression <- read.table(expression_file, sep="\t", header=TRUE)
expression <- expression[expression$phenotype_id == "ENSG00000126709",]

#jätan phenotype_id välja
exp_data <- expression[,2:length(expression)]
#sorteerin tulbad tähestiku järjekorras
exp_data <- exp_data %>% 
  dplyr::select(sort(names(.)))
exp_data <- qnorm((rank(exp_data, na.last = "keep") - 0.5) / sum(!is.na(exp_data)))
exp_data <- as.matrix(exp_data)
exp_data <- t(exp_data)
rownames(exp_data) <- c("ENSG00000126709")
colnames(exp_data) <- kept_samples_by_genotype_id$genotype_id

gene_pos <- data_frame(geneidlib="ENSG00000126709", chr="1", left=27666064, right=27672198)

covariates <- read.table(covariates_file, sep="\t", header=TRUE)
covariates_cols <- colnames(covariates)

#juhtvariandi lisamine kovariaatide hulka
juhtvariant <- read.table(juhtvariant_file, sep="\t", header=TRUE)
juhtvariant <- juhtvariant[-c(1:4)]

covariates <- rbind(covariates, juhtvariant)
covariates <- as.matrix(covariates)
colnames(covariates) <- covariates_cols
substr(colnames(covariates), 4, 4) <- "@"

runMatrixEQTL(exp_data, geno_data, snpspos, gene_pos, covariates=covariates, cisDist = 5e5, outputFileName = output_file, pvOutputThreshold = 1, permute = FALSE, model = modelLINEAR)
