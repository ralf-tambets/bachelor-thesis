library("dplyr")
library("tidyr")
library("MatrixEQTL")
source("../MatrixEQTL/MatrixEQTL_wrapper.R")

phenotype_file <- "---"
covariates_file <- "---"
#ff14_selected.R abil toodetud fail
sample_file <- "ff14_selected.tsv"
#vcfist_snpspos.txt abil toodetud fail
snpspos_file <- "ff14_snps_all.tsv"
#vcfist_doosimaatriks.txt abil toodetud fail
genotypes_file <- "ff14.dose.tsv.gz"
#ff14_expression_chunk.R abil toodetud fail
expression_file <- "ff14_expression.tsv"

kept_samples <- read.table(sample_file, sep="\t", header=TRUE)
kept_samples_by_genotype_id <- kept_samples[order(kept_samples$genotype_id),]

phenotypes <- read.table(phenotype_file, sep="\t", header = TRUE)
phenotype_gene_id <- phenotypes %>% dplyr::select(c("phenotype_id", "gene_id"))

snpspos <- read.table(snpspos_file, sep="\t", header=TRUE, row.names = NULL)
colnames(snpspos) <- c("snpid", "chr", "pos")
snpspos$chr <- as.character(snpspos$chr)

genotypes = read.table(genotypes_file, sep="\t", header=TRUE)
geno_names = paste("chr",genotypes$CHROM, "_", genotypes$POS, "_", genotypes$REF, "_", genotypes$ALT, sep="")

geno_data <- genotypes[,colnames(genotypes) %in% kept_samples$genotype_id]
#sorteerin tulbad tähestiku järjekorras
geno_data <- geno_data %>% 
  dplyr::select(sort(names(.)))
geno_data <- as.matrix(geno_data)
colnames(geno_data) <- kept_samples_by_genotype_id$sample_id
rownames(geno_data) <- geno_names

expression <- read.table(expression_file, sep="\t", header=TRUE)
expression <- left_join(expression, phenotype_gene_id, by="phenotype_id", copy=FALSE)
expression <- expression[expression$gene_id == "ENSG00000126709",]

#jätan phenotype_id ja gene_id välja
exp_data <- expression[,2:(length(expression)-1)]
#sorteerin tulbad tähestiku järjekorras
exp_data <- exp_data %>% 
  dplyr::select(sort(names(.)))
exp_data <- as.matrix(exp_data)
rownames(exp_data) <- expression$gene_id
rownames(exp_data) <- c("ENSG00000126709_a", "ENSG00000126709_b")
colnames(exp_data) <- kept_samples_by_genotype_id$sample_id

gene_pos <- data_frame(geneid=c("ENSG00000126709_a", "ENSG00000126709_b"), chr=c("1","1"), left=c(27666064,27666064), right=c(27672198,27672198))

covariates <- read.table(covariates_file, sep="\t", header=TRUE)
covariates <- covariates[,2:(length(covariates))]
colnames(covariates) <- paste("FF14", substring(colnames(covariates), first=2), sep="_")
covariates <- covariates[,colnames(covariates) %in% kept_samples$genotype_id]
covariates <- covariates %>% 
  dplyr::select(sort(names(.)))

#juhtvariandi kovariaatide hulka lisamine
juhtvariant <- read.table("juhtvariant.txt", sep="\t", header=TRUE)
juhtvariant <- juhtvariant[-c(1:4)]
juhtvariant <- juhtvariant[,colnames(juhtvariant) %in% kept_samples$genotype_id]
juhtvariant <- juhtvariant %>% 
  dplyr::select(sort(names(.)))

covariates <- rbind(covariates, juhtvariant)
covariates <- as.matrix(covariates)
colnames(covariates) <- kept_samples_by_genotype_id$sample_id

rm(kept_samples, phenotypes, genotypes, geno_names, expression)

runMatrixEQTL(exp_data, geno_data, snpspos, gene_pos, covariates=covariates, cisDist = 5e5, outputFileName = "ff14_results_IFI6_whole_genome_extra_covariate.txt", pvOutputThreshold = 1, permute = FALSE, model = modelLINEAR, noFDRsaveMemory = FALSE)
