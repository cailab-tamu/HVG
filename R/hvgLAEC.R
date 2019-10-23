library(Matrix)
library(Seurat)
library(phateR)
library(corrplot)
library(ggplot2)
library(ggrepel)

# Reading count matrix
cMatrix <- readMM('../Data/LAEC//matrix.mtx')
rownames(cMatrix) <- readLines('../Data/LAEC/genes.tsv')
colnames(cMatrix) <- readLines('../Data/LAEC/barcodes.tsv')

# Quality Control
mtRate <- colSums(cMatrix[grepl('MT-',rownames(cMatrix)),])/colSums(cMatrix)
cMatrix <- cMatrix[,mtRate < 0.1]

# Cell Cycle Assignation
cMatrix <- CreateSeuratObject(cMatrix)
cMatrix <- NormalizeData(cMatrix)
cMatrix <- ScaleData(cMatrix)
cMatrix <- CellCycleScoring(cMatrix, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
g1Cells <- cMatrix$Phase %in% 'G1'

# Filtering G1 Cells
cMatrix <- cMatrix@assays$RNA@counts
cMatrix <- cMatrix[,g1Cells]

# Normalization RC
cMatrix <- t(t(cMatrix)/colSums(cMatrix)) * 1e6

# Filtering genes with 0 counts
cMatrix <- cMatrix[rowSums(cMatrix) > 0,]

# Removing immune genes
immuneGenes <- read.csv('../Annotations/GO_term_summary_20191022_094850.txt', sep = '\t', row.names = NULL, stringsAsFactors = FALSE)
immuneGenes <- toupper(unique(immuneGenes$MGI.Gene.Marker.ID))
cMatrix <- cMatrix[!rownames(cMatrix) %in% immuneGenes,]

# Phate
drPhate <- phate(t(cMatrix), ndim = 3)
drPhate <- drPhate$embedding

# KNN of 1000 more similar cells
K <- as.matrix(dist(drPhate))
K <- t(apply(K,1,function(X){as.numeric(rank(X, 'random') <= 1000)}))

# Sampling 10 cells as center
set.seed(1)
nBoot <- 10
sCells <- sample(seq_len(ncol(cMatrix)),nBoot, replace = TRUE)

# Identify the HVG in each case
HVG <- pbapply::pbsapply(sCells, function(X){
  temp <- cMatrix[,which(as.logical(K[X,]))]
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
  #plot(xplot, col = ifelse(rownames(temp) %in% names(pAdj[pAdj < 0.01]),'red','black'))
  names(pAdj[pAdj < 0.01])
})

# Computing Jaccard Coefficient
jaccardMatrix <- sapply(HVG, function(X){
  sapply(HVG, function(Y){
    length(intersect(X,Y))/length(unique(c(X,Y)))
  })
})
newCol <- RColorBrewer::brewer.pal(9,'Blues')
png('../Results/figures/LAEC_Jaccard.png', width = 1500, height = 1500, res = 300)
corrplot::corrplot(jaccardMatrix, col = newCol, is.corr = FALSE, type = 'upper', 
                   method = 'color', mar = c(1,1,1,1), order='original',tl.col = "black")
dev.off()

# Finding union and intersection
HVG <- sort(table(unlist(HVG)), decreasing = TRUE)
HVG2 <- names(HVG[(HVG/nBoot) == 1])
outHVG <- data.frame(HVG/nBoot)
colnames(outHVG) <- c('GENE', 'FREQUENCY')
write.csv(outHVG, '../Results/hvgList/LAEC_HVG.csv')
HVG <- names(HVG)

# GO Enrichment analysis
source('https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/hsa_GO_SYMBOL.R')
GO <- hsa_GO_SYMBOL(HVG)
GO2 <- hsa_GO_SYMBOL(HVG2)
GO <- GO[order(GO$p.adjust,decreasing = FALSE),]
GO2 <- GO2[order(GO2$p.adjust,decreasing = FALSE),]
write.csv(GO, file = '../Results/annotations/LAEC_GO_UNION.csv')
write.csv(GO2, file = '../Results/annotations/LAEC_GO_INTERSECTION.csv')


# Plot
temp <- cMatrix
temp <- temp[rowMeans(temp!=0) > 0.05,]
means <- rowMeans(temp, na.rm = T)
vars <- apply(temp, 1, var, na.rm=T)
cv2 <- vars/(means^2)
recip.means <- 1/means
recip.means[is.infinite(recip.means)] <- 0
fit <- glm(cv2~recip.means, family = Gamma(link = 'identity'))
pFit <- predict(fit)
pMeans <- means[names(pFit)]
xplot <- cbind(log10(means),log10(cv2))
colnames(xplot) <- c('mean_log','cv2_log')
xplot <- as.data.frame(xplot)


source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa_ENTREZ2SYMBOL.R')
topGO <- unlist(strsplit(GO2$geneID[1], split = '/'))
topGO <- hsa_ENTREZ2SYMBOL(topGO)[,2]

png('../Results/figures/DF.png', width = 4000, height = 3000, res = 300)
gCol <- ifelse(rownames(xplot) %in% HVG, rgb(1,0,0,1), rgb(0,0,0,1))
ggplot(xplot, mapping = aes(x = mean_log, y = cv2_log)) +
  geom_point(color=gCol) +
  geom_text_repel(data = xplot[topGO,], label=topGO, max.iter = 10000, color='gray40')+
  theme_classic() +
  xlab(parse(text = 'log[10](Mean)')) +
  ylab(parse(text = 'log[10](CV^2)')) +
  ggtitle(label = 'DERMAL FIBROBLAST', subtitle = toupper(GO2$Description[1]))
dev.off()
