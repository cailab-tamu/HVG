library(Matrix)
library(phateR)
library(Seurat)
library(plot3D)
library(rgl)
library(RColorBrewer)
# Reading imputed values
igeneMatrix <- read.csv('../Data/nfkb8genes.csv')

newCol <- colorRampPalette((c('blue','skyblue','forestgreen','gold')))
par(mar=c(3,3,1,1))

sMatrix <- igeneMatrix[,c('IRF4','AICDA', 'PRDM1')]
sMatrix <- sMatrix[order(sMatrix[,3]),]
sMatrix <- as.data.frame(sMatrix)

plot(sMatrix$IRF4,sMatrix$AICDA, pch=16, cex=0.5, col=newCol(nrow(sMatrix)))
#image(t(as.matrix(sMatrix$PRDM1)), col=newCol(nrow(sMatrix)), xaxt='n')

      