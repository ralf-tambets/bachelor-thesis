library("dplyr")

expression_location <- "---"
gz.file <- gzfile(expression_location)
dat <- read.table(gz.file, sep='\t', header=TRUE)

kept_samples <- read.table("ff18_selected.tsv", sep="\t", header=TRUE)


fairfax_dat <- dat %>% 
  dplyr::select(dplyr::matches(c("phenotype_id",kept_samples$sample_id)))

write.table(fairfax_dat, file="ff18_expression.tsv", sep="\t", quote=FALSE, row.names = FALSE)