library("ggplot2")

filename <- "---"

results = read.table(filename, header=TRUE)

#histogram
ggplot(data=results, aes(p.value)) + geom_histogram()

