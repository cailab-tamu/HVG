library(Matrix)
library(phateR)
library(Seurat)
library(RColorBrewer)

# Reading imputed values
igeneMatrix <- read.csv('../Data/nfkb8genes.csv')

newCol <- colorRampPalette((c('blue','skyblue','forestgreen','gold')))
par(mar=c(3,3,1,1))

sMatrix <- igeneMatrix[,c('IRF4','AICDA', 'PRDM1')]
sMatrix <- sMatrix[order(sMatrix[,3]),]
sMatrix <- as.data.frame(sMatrix)

#Figure
png('../Results/figures/Fig4B.png', width = 2500, height = 2000, res = 300, pointsize = 20)
layout(matrix(c(1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2), nrow = 3, byrow = TRUE))
par(mar=c(3,3,1,0), mgp=c(1.5,0.5,0))
plot(sMatrix$IRF4,sMatrix$AICDA, pch=16, cex=0.6, col=newCol(nrow(sMatrix)), xlab = '',ylab='', main = 'PRDM1 (Blimp1)')
box(lwd=3)
par(mar=c(6,2.5,3,0.5), mgp=c(1.5,0.5,0))
image(t(as.matrix(sMatrix$PRDM1)), col=newCol(nrow(sMatrix)), xaxt='n', yaxt='n')
newLabels <- round(seq(min(sMatrix$PRDM1), max(sMatrix$PRDM1),(max(sMatrix$PRDM1)-min(sMatrix$PRDM1))/5),1)
axis(2,at = seq(0,1,0.2), labels = newLabels, las=2)
dev.off()      
