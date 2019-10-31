library(Matrix)
library(Seurat)
library(phateR)

# ImmuneGenes
immuneGenes <- read.csv('../Annotations/GO_term_summary_20191022_094850.txt', sep = '\t', row.names = NULL, stringsAsFactors = FALSE)
immuneGenes <- toupper(unique(immuneGenes$MGI.Gene.Marker.ID))

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

# Find HVG
findHVG <- function(X, cutOff = 0.01){
  nCells <- ncol(X)
  cellBarcodes <- colnames(X)
  phateDR <- phate(t(as.matrix(X)), ndim = 3)
  phateDR <- phateDR$embedding
  dMatrix <- as.matrix(dist(phateDR))
  set.seed(1)
  seedCell <- sample(seq_len(nCells), 1)
  seedCell <- cellBarcodes[seedCell]
  selCells <- rank(dMatrix[seedCell,]) <= 1000
  temp <- X[,selCells]
  temp <- temp[rowMeans(temp!=0) > 0.05,]
  means <- rowMeans(temp, na.rm = T)
  vars <- apply(temp, 1, var, na.rm=T)
  cv2 <- vars/(means^2)
  recip.means <- 1/means
  recip.means[is.infinite(recip.means)] <- 0
  fit <- glm(cv2~recip.means, family = Gamma(link = 'identity'))
  pFit <- predict(fit)
  pVal <- pchisq(cv2/pFit*999,999, lower.tail = FALSE)
  FC <- log2(cv2/pFit)
  pAdj <- p.adjust(pVal, method = 'fdr')
  hvgStat <- cbind(means, cv2, pFit,FC, pVal, pAdj)
  hvgStat <- as.data.frame(hvgStat)
  colnames(hvgStat) <- c('mean', 'cv2obs', 'cv2exp','log2FC', 'p.value', 'p.adj')
  HVG <- names(pAdj[pAdj < cutOff & FC > log2(1.5)])
  out <- list()
  out$dr <- phateDR
  out$seedCell <- seedCell
  out$selCells <- selCells
  out$HVG <- HVG
  out$stat <- hvgStat
  return(out)
}

# plotHVG
plotHVG <- function(X, mainLabel){
  gCol <- ifelse(rownames(X$stat) %in% X$HVG, yes = 'dodgerblue4', no = 'black')
  gPCH <- ifelse(rownames(X$stat) %in% sharedGenes, yes = 8, no = 16)
  plot(log(X$stat$mean),log(X$stat$cv2obs), col = gCol, pch = gPCH, main = mainLabel,
       xlab=parse(text = 'log(Mean)'), ylab = parse(text = 'log(CV^2)'), cex = 0.5)
  gammaReg <- X$stat[order(X$stat$mean),]
  points(log(gammaReg$mean), log(1.5*gammaReg$cv2exp), type = 'l', col = 'orange', lty=2)
  points(log(gammaReg$mean), log(gammaReg$cv2exp), type = 'l', col = 'orange')
  points(log(gammaReg$mean), log(0.66*gammaReg$cv2exp), type = 'l', col = 'orange', lty=2)
  
  HVG <- length(X$HVG)
  noHVG <- nrow(X$stat)-HVG
  nShared <- length(sharedGenes)
  legend('topright', legend = c(paste0('No HVG (',noHVG,')'), paste0('HVG (',HVG,')'), paste0('Shared HVG (',nShared,')'), 'Gamma\nRegression'), col = c('black', 'dodgerblue4', 'dodgerblue4', 'orange'),pch = c(16,16,8, NA), lty = c(NA,NA,NA,1),  bty = 'n')
}

# v2_hIPSC
v2IPSC <- readMM('../Data/v2_hIPSC/SRS2119548.mtx')
colnames(v2IPSC) <- readLines('../Data/v2_hIPSC/barcodes.tsv')
rownames(v2IPSC) <- readLines('../Data/v2_hIPSC/genes.tsv')
v2IPSC <- scQC(v2IPSC)
hvgV2 <- findHVG(v2IPSC)

# s4_hIPSC
s4IPSC <- readMM('../Data/s4_hIPSC/E-MTAB-6687.mtx')
rownames(s4IPSC) <- readLines('../Data/s4_hIPSC/genes.tsv')
colnames(s4IPSC) <- readLines('../Data/s4_hIPSC/barcodes.tsv')
s4IPSC <- scQC(s4IPSC)
hvgS4 <- findHVG(s4IPSC)

# Plot
png(filename = '../Results/figures/allHIPSC.png', width = 4000, height = 2000, res = 300, pointsize = 15)
sharedGenes <- intersect(hvgS4$HVG, hvgV2$HVG)
par(mfrow=c(1,2), mar=c(3,3,2,1), mgp=c(1.5,0.5,0))
plotHVG(hvgV2, mainLabel = 'hiPSC V2 INTERNEURONS')
plotHVG(hvgS4, mainLabel = 'hiPSC FIBROBLASTS')
dev.off()
