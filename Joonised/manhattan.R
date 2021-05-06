library(qqman)
library(reshape2)

metaanalysis_file <- "---"

meta_results = read.table(metaanalysis_file, header=TRUE)

osad = colsplit(string=meta_results$MarkerName, pattern="_", names=c("chr", "bp", "ref", "alt"))

osad$chr <- gsub("chr", "", osad$chr)

manhattan_df <- data.frame(SNP=meta_results$MarkerName, CHR=as.numeric(osad$chr), BP=osad$bp, P=meta_results$P.value)

rm(meta_results)

rm(osad)

manhattan(manhattan_df, main="Manhattan plot, 4 studies", cex=0.8, cex.axis=0.9, col=c("blue4", "orange3"))
