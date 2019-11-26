library(Matrix)
library(Seurat)
library(phateR)

# ImmuneGenes
immuneGenes <- read.csv('../Annotations/GO_term_summary_20191022_094850.txt', sep = '\t', row.names = NULL, stringsAsFactors = FALSE)
immuneGenes <- toupper(unique(immuneGenes$MGI.Gene.Marker.ID))

# Filter by library size, mitochondrial rate and cell-cycle phase and CPM normalization
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
  pVal <- pchisq((cv2/pFit)*999,999, lower.tail = FALSE)
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

# Plot
plotHVG <- function(X, mainLabel, ...){
  gCol <- ifelse(rownames(X$stat) %in% X$HVG, yes = 'dodgerblue4', no = 'black')
  gPCH <- 20 #ifelse(rownames(X$stat) %in% sharedGenes, yes = 8, no = 16)
  plot(log(X$stat$mean),log(X$stat$cv2obs), col = gCol, pch = gPCH, main = mainLabel,
       xlab=parse(text = 'log(Mean)'), ylab = parse(text = 'log(CV^2)'), cex = 0.5,...)
  gammaReg <- X$stat[order(X$stat$mean),]
  points(log(gammaReg$mean), log(1.5*gammaReg$cv2exp), type = 'l', col = 'orange', lty=2)
  points(log(gammaReg$mean), log(gammaReg$cv2exp), type = 'l', col = 'orange')
  points(log(gammaReg$mean), log(0.66*gammaReg$cv2exp), type = 'l', col = 'orange', lty=2)
  
  HVG <- length(X$HVG)
  noHVG <- nrow(X$stat)-HVG
  legend('topright', legend = c('HVG', 'Gamma\nRegression'), col = c('dodgerblue4', 'orange'),pch = c(16, NA), lty = c(NA,1),  bty = 'n', pt.cex = 0.5)
}

# iPSC
iPSC <- readMM('../Data/s4_hiPSC/E-MTAB-6687.mtx')
colnames(iPSC) <- readLines('../Data/s4_hiPSC/barcodes.tsv')
rownames(iPSC) <- readLines('../Data/s4_hiPSC/genes.tsv')
iPSC <- scQC(iPSC)
iPSC <- findHVG(iPSC)

# B-Cells
B <- readMM('../Data/GM12878/matrix.mtx')
rownames(B) <- read.table('../Data/GM12878/genes.tsv', stringsAsFactors = FALSE)[,2]
colnames(B) <- readLines('../Data/GM12878/barcodes.tsv')
B <- scQC(B)
B <- findHVG(B)

# LAEC
LAEC <- readMM('../Data/LAEC/matrix.mtx')
rownames(LAEC) <- readLines('../Data/LAEC/genes.tsv')
colnames(LAEC) <- readLines('../Data/LAEC/barcodes.tsv')
LAEC <- scQC(LAEC)
LAEC <- findHVG(LAEC)

# DF
DF <- readMM('../Data/DF/matrix.mtx')
rownames(DF) <- read.table('../Data/DF/genes.tsv', stringsAsFactors = FALSE)[,2]
colnames(DF) <- readLines('../Data/DF/barcodes.tsv')
DF <- scQC(DF)
DF <- findHVG(DF)

# Figure 6
png('../Results/figures/Fig6.png', width = 6000, height = 1500, res = 300, pointsize = 20)
par(mfrow=c(1,4), mar=c(3,3,2.5,1), mgp=c(1.5,0.5,0))
plotHVG(iPSC, 'iPSC')
plotHVG(B, 'LCL')
plotHVG(LAEC, 'LUNG AIRWAY EPITHELIUM')
plotHVG(DF, 'DERMAL FIBROBLAST')
dev.off()
