library(Matrix)
library(phateR)
library(Seurat)
library(plot3D)
library(rgl)

# Phate DR
drPhate <- read.csv('../Data/GM12878/GM12878PhateCoordinates.csv', header = FALSE)
drPhate <- drPhate[order(drPhate$V1, decreasing = TRUE),]

# Reading count matrix
GM12878 <- readMM('../Data/GM12878/matrix.mtx')
rownames(GM12878) <- read.table('../Data/GM12878/genes.tsv', row.names = 1, stringsAsFactors = FALSE)[,1]
colnames(GM12878) <- readLines('../Data/GM12878/barcodes.tsv')

# Phate DR
GM12878 <- t(t(GM12878)/colSums(GM12878))*1e6
drPhate2 <- phate(as.matrix(t(GM12878)), ndim = 3)
drPhate2 <- drPhate2$embedding

png('../Results/figures/FIG1.png', width = 6000, height = 2000, res = 300, pointsize = 20)

# Plot A
lDistribution <- matrix(c(1,1,1,2,3,4,5,6,1,1,1,7,8,9,10,11), nrow = 2, byrow = TRUE)
layout(lDistribution)
par(mar=c(4,3,1,0), mgp=c(1.5,0.5,0))

cellColor <- densCols(drPhate)

scatter3D(drPhate[,1],drPhate[,2], drPhate[,3], pch = 16, theta = 0, bty = 'b',
          phi = 0, colvar = FALSE, cex = 0.5, col = cellColor, xlab = 'PHATE 1', cex=0.5,
          ylab = 'PHATE 2', zlab = 'PHATE 3', main = 'LCL - GM12878',ticktype = "detailed")
scatter3D(x = -10,y = -100, z = 0, cex=14, add = TRUE, col='red', lwd = 5)

# Plot B
K <- as.matrix(dist(drPhate2))
sCells <- K[rownames(drPhate2) %in% 'ACGATGTTCTAACTCT-1',]
sCells <- rank(sCells) <= 1000
sCells <- GM12878[,sCells]
sCells <- CreateSeuratObject(sCells)
sCells <- NormalizeData(sCells)
sCells <- ScaleData(sCells)
sCells <- FindVariableFeatures(sCells)
sCells <- RunPCA(sCells, verbose = FALSE)
for(i in seq(10,100,10)){
  sCells <- RunTSNE(sCells, perplexity = i, check_duplicates=FALSE)
  tPositions <- sCells@reductions$tsne@cell.embeddings
  par(mar=c(2.5,2.5,2,1), mgp=c(1.5,0.5,0))
  plot(tPositions, xlab = 't-SNE 1', ylab = 't-SNE 2', col = densCols(tPositions), pch = 16, cex = 0.5, main = parse(text = paste0('Perplexity:',i)))
}
dev.off()
