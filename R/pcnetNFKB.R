library(Matrix)
library(Seurat)

GM12878 <- readMM('../Data/GM12878/matrix.mtx')
rownames(GM12878) <- read.table('../Data/GM12878/genes.tsv', stringsAsFactors = FALSE)[,2]
colnames(GM12878) <- readLines('../Data/GM12878/barcodes.tsv')

scQC <- function(X){
  X <- X[,colSums(X) > 1000]
  lSize <- colSums(X)
  X <- X[,!lSize %in% boxplot.stats(lSize)$out]
  mtRate <- colSums(X[grepl('MT-',rownames(X)),])/colSums(X)
  X <- X[,mtRate < 0.1]
  X <- CreateSeuratObject(X)
  X <- NormalizeData(X)
  X <- ScaleData(X)
  X <- CellCycleScoring(X, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  X <- X@assays$RNA@counts[,X$Phase == 'G1']
  X <- X[rowSums(X) > 0,]
  return(X)
}

GM12878 <- scQC(GM12878)


eAssociations <- read.csv('../Data/nfkbAssociations.csv', row.names = 1)

gExpression <- GM12878[rownames(eAssociations),]

source('https://raw.githubusercontent.com/cailab-tamu/scTenifoldNet/master/R/pcNet.R')
oAssociation <- pcNet(gExpression, nComp = 3, symmetric = TRUE, scaleScores = TRUE) 
oAssociation[is.na(eAssociations)] <- NA
oAssociation <- as.matrix(oAssociation)

library(corrplot)
corrplot.mixed(oAssociation, lower.col = 'black', na.label = ' ', mar = c(1,1,1,1))
  