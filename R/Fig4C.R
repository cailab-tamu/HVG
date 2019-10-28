library(Matrix)
library(phateR)
library(Seurat)
library(plot3D)
library(rgl)

# Reading count matrix
GM12878 <- readMM('../Data/GM12878/matrix.mtx')
rownames(GM12878) <- read.table('../Data/GM12878/genes.tsv', row.names = 1, stringsAsFactors = FALSE)[,1]
colnames(GM12878) <- readLines('../Data/GM12878/barcodes.tsv')

# Imputation
library(Rmagic)
iGM12878 <- t(magic(t(as.matrix(GM12878)))$result)

# QC Library size and mitochondrial rate
lSize <- colSums(GM12878)
GM12878 <- GM12878[,names(lSize[lSize > 1000])]
lSize <- colSums(GM12878)
GM12878 <- GM12878[,names(lSize[!lSize %in% boxplot.stats(lSize)$out])]
mtRate <- colSums(GM12878[grepl('MT-', rownames(GM12878)),])/colSums(GM12878)
GM12878 <- GM12878[,names(mtRate[mtRate < 0.1])]

# Cell cycle assignation
GM12878 <- CreateSeuratObject(GM12878)
GM12878 <- NormalizeData(GM12878)
GM12878 <- ScaleData(GM12878)
GM12878 <- CellCycleScoring(GM12878,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
g1Cells <- GM12878$Phase %in% 'G1'

# Filtering G1 Cells
GM12878 <- GM12878@assays$RNA@counts
GM12878 <- GM12878[,g1Cells]
iGM12878 <- iGM12878[,g1Cells]

#Normalization
GM12878 <- t(t(GM12878)/colSums(GM12878))*1e6

# Genes
geneList <- c('AICDA', 'BACH2', 'BCL6', 'IRF4', 'PAX5', 'PRDM1', 'REL','RELA')


geneMatrix <- GM12878[geneList,]
igeneMatrix <- iGM12878[geneList,]
cor(as.matrix(t(igeneMatrix)), method = 'sp')

png('../Results/figures/FIG4C.png', width = 2000, height = 2000, res = 300, pointsize = 20)
par(mfrow=c(8,8), mar=c(.5,.5,.5,.5))
for(i in 1:8){
  for(j in 1:8){
    if(i == j){
        hist((igeneMatrix[i,]), main = '',xaxt='n', yaxt='n')
      if(i == 1 & j == 1){
        mtext(text = geneList[j], line = -0.1, side = 3, cex = 0.4)
        mtext(text = geneList[i], line = -0.1, side = 2, cex = 0.4)
      }
    } else {
      if(i == 1){
        data <- cbind(igeneMatrix[i,], igeneMatrix[j,])
        data <- densCols(data)
        plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = data)
        abline(lm(igeneMatrix[j,]~igeneMatrix[i,]), col='red')
        mtext(text = geneList[j], line = -0.1, side = 3, cex = 0.4)
        corVal <- cor(igeneMatrix[i,], igeneMatrix[j,], method = 'spearman')
        legend('topleft', legend = parse(text = paste0('rho ==', round(corVal,2))), bty = 'n', cex = 0.4)
        
      }else{
        if(j == 1){
          data <- cbind(igeneMatrix[i,], igeneMatrix[j,])
          data <- densCols(data)
          plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = data)
          abline(lm(igeneMatrix[j,]~igeneMatrix[i,]), col='red')
          mtext(text = geneList[i], line = -0.1, side = 2, cex = 0.4)
          corVal <- cor(igeneMatrix[i,], igeneMatrix[j,], method = 'spearman')
          legend('topleft', legend = parse(text = paste0('rho ==', round(corVal,2))), bty = 'n', cex = 0.4)
        }else {
          data <- cbind(igeneMatrix[i,], igeneMatrix[j,])
          data <- densCols(data)
          plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = data)
          abline(lm(igeneMatrix[j,]~igeneMatrix[i,]), col='red')
          corVal <- cor(igeneMatrix[i,], igeneMatrix[j,], method = 'spearman')
          legend('topleft', legend = parse(text = paste0('rho ==', round(corVal,2))), bty = 'n', cex = 0.4)
        }
      }
    }
  }
}
dev.off()

