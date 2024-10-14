## start R studio session
setwd("C:/data/sampleQC_preQC/FPCA/run1")
list.files(path="./")   ##list files in current directory
library(RColorBrewer)

##*sa file tab delimited with headers: sample-id       sex     affection       family-id       father-id       mother-id
##*sa file is produced after running fpca step (FRATOOLS package)
sample_annotation <- read.table("samples.sa", header=T, sep = "\t")
data <- read.table("sieved-exome_samples_removed.common-pruned-recode_qced.pca", header=T, sep = "\t")
data <- merge(sample_annotation, data, by = "sample.id")


### PC1 vs PC2 ###
x11(10,7, pointsize = 13)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
col = brewer.pal(5,"RdBu")
col = c("red", "yellow", "blue", "white", "cyan", "purple", "green","pink","black","cyan","magenta","orange")
#col = c(col[1], col[3], "green","purple", col[5])

plot(data$PC1, data$PC2, xlab="PC1", ylab="PC2", main="PC1 vs PC2", pch = 21, bg = col[unclass(factor(data$affection))])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$affection))), pt.bg = col,
        text.col = "black", pch = 21,
        bg = 'white')


### PC3 vs PC4 ###
x11(10,7, pointsize = 13)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
col = brewer.pal(5,"RdBu")
col = c("red", "yellow", "blue", "white", "cyan", "purple", "green","pink","black","cyan","magenta","orange")
#col = c(col[1], col[3], "green","purple", col[5])

plot(data$PC3, data$PC4, xlab="PC3", ylab="PC4", main="PC3 vs PC4", pch = 21, bg = col[unclass(factor(data$affection))])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$affection))), pt.bg = col, text.col = "black", pch = 21, bg = 'white')



### PC5 vs PC6 ###
x11(10,7, pointsize = 13)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
col = brewer.pal(5,"RdBu")
col = c("red", "yellow", "blue", "white", "cyan", "purple", "green","pink","red","black","magenta","orange")
#col = c(col[1], col[3], "green","purple", col[5])

plot(data$PC5, data$PC6, xlab="PC5", ylab="PC6", main="PC5 vs PC6", pch = 21, bg = col[unclass(factor(data$affection))])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$affection))), pt.bg = col, text.col = "black", pch = 21, bg = 'white')

