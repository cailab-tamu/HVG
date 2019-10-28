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
  pAdj <- p.adjust(pVal, method = 'fdr')
  hvgStat <- cbind(means, cv2, pFit, pVal, pAdj)
  hvgStat <- as.data.frame(hvgStat)
  colnames(hvgStat) <- c('mean', 'cv2obs', 'cv2exp', 'p.value', 'p.adj')
  HVG <- names(pAdj[pAdj < cutOff])
  out <- list()
  out$dr <- phateDR
  out$seedCell <- seedCell
  out$selCells <- selCells
  out$HVG <- HVG
  out$stat <- hvgStat
  return(out)
}

# DERMAL FIBROBLASTS
DF <- readMM('../Data/DF/matrix.mtx')
rownames(DF) <- read.csv('../Data/DF/genes.tsv', sep = '\t', stringsAsFactors = FALSE, header = FALSE)[,2]
colnames(DF) <- readLines('../Data/DF/barcodes.tsv')
DF <- scQC(DF)
hvgDF <- findHVG(DF)

# LUNG DISTAL FIBROBLASTS
LDF <- readMM('../Data/FIBROBLAST1/SRS2769050.mtx')
rownames(LDF) <- readLines('../Data/FIBROBLAST1/genes.tsv')
colnames(LDF) <- readLines('../Data/FIBROBLAST1/barcodes.tsv')
LDF <- scQC(LDF)
hvgLDF <- findHVG(LDF)

# LUNG PROXIMAL FIBROBLASTS
LPF <- readMM('../Data/FIBROBLAST2/SRS2769051.mtx')
rownames(LPF) <- readLines('../Data/FIBROBLAST2/genes.tsv')
colnames(LPF) <- readLines('../Data/FIBROBLAST2/barcodes.tsv')
LPF <- scQC(LPF)
hvgLPF <- findHVG(LPF)

png('../Results/figures/allFibroblasts.png', width = 6000, height = 2000, res = 300, pointsize = 20)
par(mfrow=c(1,3), mar=c(2.5,3,2,1), mgp=c(1.5,0.5,0))
sharedGenes <- intersect(intersect(hvgDF$HVG, hvgLDF$HVG), hvgLPF$HVG)
intersect(hvgLDF$HVG, hvgLPF$HVG)
gCol <- ifelse(hvgDF$stat$p.adj < 0.01, yes = 'dodgerblue4', no = 'black')
gPCH <- ifelse(rownames(hvgDF$stat) %in% sharedGenes, yes = 8, no = 16)
plot(log(hvgDF$stat$mean),log(hvgDF$stat$cv2obs), col = gCol, pch = gPCH, main = 'DERMAL FIBROBLAST',
     xlab=parse(text = 'log(Mean)'), ylab = parse(text = 'log(CV^2)'))
gammaReg <- hvgDF$stat[order(hvgDF$stat$mean),]
points(log(gammaReg$mean), log(gammaReg$cv2exp), type = 'l', col = 'orange')
legend('topright', legend = c('No HVG', 'HVG FDR < 0.01', 'Shared HVG', 'Gamma\nRegression'), col = c('black', 'dodgerblue4', 'dodgerblue4', 'orange'),pch = c(16,16,8, NA), lty = c(NA,NA,NA,1),  bty = 'n')

gCol <- ifelse(hvgLDF$stat$p.adj < 0.01, yes = 'dodgerblue4', no = 'black')
gPCH <- ifelse(rownames(hvgLDF$stat) %in% sharedGenes, yes = 8, no = 16)
plot(log(hvgLDF$stat$mean),log(hvgLDF$stat$cv2obs), col = gCol, pch = gPCH, main = 'LUNG DISTAL FIBROBLAST',
     xlab=parse(text = 'log(Mean)'), ylab = parse(text = 'log(CV^2)'))
gCol <- ifelse(hvgLPF$stat$p.adj < 0.01, yes = 'dodgerblue4', no = 'black')
gammaReg <- hvgLDF$stat[order(hvgLDF$stat$mean),]
points(log(gammaReg$mean), log(gammaReg$cv2exp), type = 'l', col = 'orange')
legend('topright', legend = c('No HVG', 'HVG FDR < 0.01', 'Shared HVG', 'Gamma\nRegression'), col = c('black', 'dodgerblue4', 'dodgerblue4', 'orange'),pch = c(16,16,8, NA), lty = c(NA,NA,NA,1),  bty = 'n')

gPCH <- ifelse(rownames(hvgLPF$stat) %in% sharedGenes, yes = 8, no = 16)
plot(log(hvgLPF$stat$mean),log(hvgLPF$stat$cv2obs), col = gCol, pch = gPCH, main = 'LUNG PROXIMAL FIBROBLAST',
     xlab=parse(text = 'log(Mean)'), ylab = parse(text = 'log(CV^2)'))
gammaReg <- hvgLPF$stat[order(hvgLPF$stat$mean),]
points(log(gammaReg$mean), log(gammaReg$cv2exp), type = 'l', col = 'orange')
legend('topright', legend = c('No HVG', 'HVG FDR < 0.01', 'Shared HVG', 'Gamma\nRegression'), col = c('black', 'dodgerblue4', 'dodgerblue4', 'orange'),pch = c(16,16,8, NA), lty = c(NA,NA,NA,1),  bty = 'n')
dev.off()

source('https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/hsa_GO_SYMBOL.R')
GO <- hsa_GO_SYMBOL(sharedGenes)
GO <- GO[order(GO$p.adjust),]
write.csv(GO, '../Results/annotations/allFIBROBLASTS.csv')
