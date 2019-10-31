library(Matrix)
library(Seurat)
library(phateR)
library(corrplot)
library(ggplot2)
library(ggrepel)

# Reading count matrix
cMatrix <- readMM('../Data/DF/matrix.mtx')
rownames(cMatrix) <- read.csv('../Data/DF/genes.tsv', header = FALSE, sep = '\t', stringsAsFactors = FALSE)[,2]
colnames(cMatrix) <- readLines('../Data/DF/barcodes.tsv')

# Loading immune genes
immuneGenes <- read.csv('../Annotations/GO_term_summary_20191022_094850.txt', sep = '\t', row.names = NULL, stringsAsFactors = FALSE)
immuneGenes <- toupper(unique(immuneGenes$MGI.Gene.Marker.ID))

# Quality Control
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
  X <- (t(t(X)/colSums(X))*1e6)
  X <- X[rowSums(X) > 0,]
  X <- X[!rownames(X) %in% immuneGenes,]
  X <- X[!grepl('MT-',rownames(X)),]
  return(X)
}
cMatrix <- scQC(cMatrix)

# Phate
drPhate <- phate(t(cMatrix), ndim = 3)
drPhate <- drPhate$embedding

# Sampling 10 cells as center
set.seed(1)
nBoot <- 10
sCells <- sample(seq_len(ncol(cMatrix)),nBoot, replace = TRUE)

# KNN of 1000 more similar cells
K <- as.matrix(dist(drPhate))
K <- K[sCells,]
K <- t(apply(K,1,function(X){rank(X) <= 1000}))

png('../Results/figures/rep10DF.png', width = 6000, height = 2000, pointsize = 15, res = 300)
lMatrix <- matrix(data = c(1,2,3,4,5,11,11,6,7,8,9,10,11,11), nrow = 2, byrow = TRUE)
layout(lMatrix)
par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))

# Identify the HVG in each case
HVG <- pbapply::pbsapply(seq_len(nBoot), function(X){
  temp <- cMatrix[,which(as.logical(K[X,]))]
  temp <- temp[rowMeans(temp!=0) > 0.05,]
  means <- rowMeans(temp, na.rm = T)
  vars <- apply(temp, 1, var, na.rm=T)
  cv2 <- vars/(means^2)

  recip.means <- 1/means
  recip.means[is.infinite(recip.means)] <- 0

  fit <- glm(cv2~recip.means, family = Gamma(link = 'identity'))
  pFit <- predict(fit)
  FC <- cv2/pFit
  pVal <- pchisq((cv2/pFit)*999,999, lower.tail = FALSE)
  pAdj <- p.adjust(pVal, method = 'fdr')
  #plot(xplot, col = ifelse(rownames(temp) %in% names(pAdj[pAdj < 0.01]),'red','black'))
  listHVG <- names(pAdj[pAdj < 0.01 & FC > 1.5])
  gCol <- ifelse(rownames(temp) %in% listHVG, 'dodgerblue4','black')
  regLine <- cbind(log(means),(pFit))
  regLine <- regLine[order(regLine[,1]),]
  plot(log(means),log(cv2), col = gCol, xlab=parse(text = 'log(Mean)'), ylab=parse(text = 'log(CV^2)'), pch = 16, cex = 0.5, main = paste0('REPLICATE: ', X))
  points(regLine[,1], log(1.5*regLine[,2]), col='orange', type = 'l', lty = 2)
  points(regLine[,1],log(regLine[,2]), col='orange', type = 'l')
  points(regLine[,1], log(0.66*regLine[,2]), col='orange', type = 'l', lty = 2)  
  nHVG <- length(listHVG)
  legend('bottomleft', legend = c(paste0('HVG (',nHVG,')'), 'Gamma\nRegression'), bty = 'n', pch = c(16,NA), col=c('dodgerblue4','orange'), lty=c(NA,1))
  return(listHVG)
})

# Computing Jaccard Coefficient
jaccardMatrix <- sapply(HVG, function(X){
  sapply(HVG, function(Y){
    length(intersect(X,Y))/length(unique(c(X,Y)))
  })
})
newCol <- RColorBrewer::brewer.pal(9,'Blues')
#png('../Results/figures/DF_Jaccard.png', width = 1500, height = 1500, res = 300)
corrplot::corrplot(jaccardMatrix, col = newCol, is.corr = FALSE, type = 'upper', 
                   method = 'color', mar = c(1,1,1,1), order='original',tl.col = "black", title = toupper('Jaccard Similarity Index'))
dev.off()

# Finding union and intersection
HVG <- sort(table(unlist(HVG)), decreasing = TRUE)
HVG2 <- names(HVG[(HVG/nBoot) == 1])
outHVG <- data.frame(HVG/nBoot)
colnames(outHVG) <- c('GENE', 'FREQUENCY')
write.csv(outHVG, '../Results/hvgList/DF_HVG.csv')
HVG <- names(HVG)

# GO Enrichment analysis
source('https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/hsa_GO_SYMBOL.R')
GO <- hsa_GO_SYMBOL(HVG)
GO2 <- hsa_GO_SYMBOL(HVG2)
GO <- GO[order(GO$p.adjust,decreasing = FALSE),]
GO2 <- GO2[order(GO2$p.adjust,decreasing = FALSE),]
write.csv(GO, file = '../Results/annotations/DF_GO_UNION.csv')
write.csv(GO2, file = '../Results/annotations/DF_GO_INTERSECTION.csv')
