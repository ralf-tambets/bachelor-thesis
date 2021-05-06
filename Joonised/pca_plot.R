library("ggplot2")

#Quach_2016 põhjal arvutatud PCA fail
filename <- "---"

pca <- read.table(filename, sep="\t", header=TRUE)

pca <- data.frame(t(pca))

colnames(pca) <- c("PC1","PC2","PC3","PC4","PC5","PC6")

ggplot(pca, aes(x=PC1, y=PC2)) + geom_point()
