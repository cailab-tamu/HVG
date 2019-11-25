library(Matrix)
library(Seurat)
library(phateR)

EMTAB6268 <- readMM('../Data/E-MTAB-6268/E-MTAB-6268Day0Rep1.mtx')
rownames(EMTAB6268) <- readLines('../Data/E-MTAB-6268/genes_E-MTAB-6268Day0Rep1.tsv')
colnames(EMTAB6268) <- readLines('../Data/E-MTAB-6268/barcodes_E-MTAB-6268Day0Rep1.tsv')

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
  set.seed(1)
  phateDR <- phate(t(as.matrix(X)), ndim = 3)
  phateDR <- phateDR$embedding
  dMatrix <- as.matrix(dist(phateDR))
  set.seed(2)
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
  length(HVG)
  out <- list()
  out$dr <- phateDR
  out$seedCell <- seedCell
  out$selCells <- selCells
  out$HVG <- HVG
  out$stat <- hvgStat
  return(out)
}

plotHVG <- function(X, mainLabel){
  gCol <- ifelse(rownames(X$stat) %in% X$HVG, yes = 'dodgerblue4', no = 'black')
  gPCH <- 16
  plot(log(X$stat$mean),log(X$stat$cv2obs), col = gCol, pch = gPCH, main = mainLabel,
       xlab=parse(text = 'log(Mean)'), ylab = parse(text = 'log(CV^2)'), cex = 0.5)
  gammaReg <- X$stat[order(X$stat$mean),]
  points(log(gammaReg$mean), log(1.5*gammaReg$cv2exp), type = 'l', col = 'orange', lty=2)
  points(log(gammaReg$mean), log(gammaReg$cv2exp), type = 'l', col = 'orange')
  points(log(gammaReg$mean), log(0.66*gammaReg$cv2exp), type = 'l', col = 'orange', lty=2)
  
  HVG <- length(X$HVG)
  noHVG <- nrow(X$stat)-HVG
  legend('topright', legend = c(paste0('No HVG (',noHVG,')'), paste0('HVG (',HVG,')'), 'Gamma\nRegression'), col = c('black', 'dodgerblue4', 'orange'),pch = c(16,16,NA), lty = c(NA,NA,1),  bty = 'n')
}

EMTAB6268 <- scQC(EMTAB6268)
hvgEMTAB6268 <- findHVG(EMTAB6268)

plotHVG(hvgEMTAB6268, 'E-MTAB-6268')
writeLines(hvgEMTAB6268$HVG, sep=', ')
# ID3, LEFTY1, RPP14, HIST1H4C, HSP90AB1, MALAT1, TAGLN, KPNA2, UBE2C, PARD6B, ZNF91, MT-ND1, MT-ND2, MT-CO1, MT-CO2, MT-ATP6, MT-CO3, MT-ND4, MT-CYB
