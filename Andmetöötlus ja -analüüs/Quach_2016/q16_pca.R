library(SNPRelate)

output_file <- "---"
#q16_vcfist_gds.R abil toodetud fail
gds_file <- "---"
genofile <- snpgdsOpen(gds_file)

#LD-based SNP pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(unname(snpset))

#PCA analysis
pca <- snpgdsPCA(genofile, snp.id=snpset.id, eigen.cnt = 6)

#pane fail kinni
snpgdsClose(genofile)

genotype.id <- pca$sample.id

tab <- data.frame(geno_PC1 = pca$eigenvect[,1],
                  geno_PC2 = pca$eigenvect[,2],
                  geno_PC3 = pca$eigenvect[,3],
                  geno_PC4 = pca$eigenvect[,4],
                  geno_PC5 = pca$eigenvect[,5],
                  geno_PC6 = pca$eigenvect[,6],
                  
                  stringsAsFactors = FALSE)
mtab <- t(tab)
colnames(mtab) <- genotype.id
rownames(mtab) <- c()

write.table(mtab, file=output_file, sep="\t", quote = FALSE, row.names = FALSE)
