library("qqman")
library("ggplot2")

filename <- "---"

results = read.table(filename, header=TRUE)

#qq-plot
qq(results$p.value)

