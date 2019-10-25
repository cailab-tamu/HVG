library(Matrix)
library(phateR)
library(Seurat)
library(plot3D)
library(rgl)

# Reading count matrix
GM12878 <- readMM('../Data/GM12878/matrix.mtx')
rownames(GM12878) <- read.table('../Data/GM12878/genes.tsv', row.names = 1, stringsAsFactors = FALSE)[,1]
colnames(GM12878) <- readLines('../Data/GM12878/barcodes.tsv')

# Phate DR
drPhate <- phate(t(GM12878), ndim = 3)
drPhate <- drPhate$embedding

# Plot A
png('../Results/figures/FIG1.png', width = 6000, height = 2000, res = 300, pointsize = 20)
lDistribution <- matrix(c(1,1,1,2,3,4,5,6,1,1,1,7,8,9,10,11), nrow = 2, byrow = TRUE)
layout(lDistribution)
par(mar=c(0,0,2,0), mgp=c(1.5,0.5,0))
K <- as.matrix(dist(drPhate))

cellColor <- densCols(drPhate)
sCells <- K[rownames(drPhate) %in% 'CCCTCCTAGGGAAACA-1',]
cellColor[rank(sCells) <= 1000] <- rgb(1,0,0,1)

scatter3D(drPhate[,1],drPhate[,2], drPhate[,3], pch = 16, theta = -45, 
          phi = 45, colvar = FALSE, col = cellColor, xlab = 'PHATE 1', cex=0.5,
          ylab = 'PHATE 2', zlab = 'PHATE 3', main = 'LCL - GM12878\nCCCTCCTAGGGAAACA')

# Plot B
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
