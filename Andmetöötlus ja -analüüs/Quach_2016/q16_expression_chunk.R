library("dplyr")

output_file <- "---"
samples_file <- "---"
expression_location <- "---"
gz.file <- gzfile(expression_location)
dat <- read.table(gz.file, sep='\t', header=TRUE)

kept_samples <- read.table(samples_file, sep="\t", header=TRUE)


fairfax_dat <- dat %>% 
  dplyr::select(dplyr::matches(c("phenotype_id",kept_samples$sample_id)))

write.table(fairfax_dat, file=output_file, sep="\t", quote=FALSE, row.names = FALSE)