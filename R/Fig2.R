library(Matrix)
library(phateR)
library(Seurat)
# library(plot3D)
# library(rgl)

# Reading count matrix
GM12878 <- readMM('../Data/GM12878/matrix.mtx')
rownames(GM12878) <- read.table('../Data/GM12878/genes.tsv', row.names = 1, stringsAsFactors = FALSE)[,1]
colnames(GM12878) <- readLines('../Data/GM12878/barcodes.tsv')

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

#Normalization
GM12878 <- t(t(GM12878)/colSums(GM12878))*1e6

# Regression model
temp <- GM12878
temp <- temp[rowMeans(temp!=0) > 0.01,]
means <- rowMeans(temp, na.rm = T)
vars <- apply(temp, 1, var, na.rm=T)
cv2 <- vars/(means^2)
recip.means <- 1/means
recip.means[is.infinite(recip.means)] <- 0
fit <- glm(cv2~recip.means, family = Gamma(link = 'identity'))
pFit <- predict(fit)
pMeans <- means[names(pFit)]
pVal <- pchisq((cv2/pFit)*999,999, lower.tail = FALSE)
pAdj <- p.adjust(pVal, 'fdr')

# Plot
png('../Results/figures/FIG2.png',width = 6000, height = 2000, res = 300, pointsize = 20)
lMatrix <- matrix(c(1,1,2,3,4,1,1,5,6,7), nrow = 2, byrow = TRUE)
layout(lMatrix)
xplot <- cbind(log(means),log(cv2))
colnames(xplot) <- c('mean_log','cv2_log')
xplot <- as.data.frame(xplot)
par(mar=c(3,3,1,1), mgp = c(1.5,0.5,0))
FC <- cv2/pFit
geneColor <- ifelse(pAdj < 0.01 & FC > 1.5, yes = 'dodgerblue4', no = 'black')
genePch <- rep(16, nrow(xplot))
genePch[rownames(xplot) %in% c('TMEM9B', 'IGKC', 'LTB', 'RPL17', 'CCL3', 'FTL')] <- 8
plot(xplot, pch = genePch, col=geneColor, 
     xlab = 'log(Mean Expression)', 
     ylab=parse(text = 'log(CV^2)'), 
     cex = 0.5)
regLine <- cbind(log(pMeans),pFit)
regLine <- regLine[order(regLine[,1]),]
points(regLine[,1], log(regLine[,2]), type='l', col= 'orange')
points(regLine[,1], log(1.5*regLine[,2]), type='l', col= 'orange', lty=2)
points(regLine[,1], log(0.66*regLine[,2]), type='l', col= 'orange', lty=2)

text(4.165297-1,0.1695131, 'TMEM9B', col = 'black')
text(8.580823+0.5,2.831446, 'IGKC', col = 'black')
text(5.542169+0.5, 2.288525, 'LTB', col = 'black')
text(6.17042-1,-1.501624, 'RPL17', col = 'black')
text(5.09993+0.5, 3.09279, 'CCL3', col = 'black')
text(8.352254+0.5,0.04420679, 'FTL', col = 'black')
legend('topright', legend = c('HVG', 'Gamma\nRegression'), bty = 'n', pch = c(16,NA), lty = c(NA,1), pt.cex = 0.5, col = c('dodgerblue4', 'orange'))


par(mar=c(2,3,1,1), mgp = c(1.5,0.5,0))
barplot(as.numeric(GM12878['TMEM9B',]), main = 'TMEM9B', ylab = 'Expression Level (CPM)', yaxt = 'n', ylim=c(0,550))
axis(2, at = seq(0, 500, 500/3), labels = formatC(seq(0, 500, 500/3), digits = 0, format = 'e'))
mtext('Cell Index', side = 1, cex = 0.7)
box()
barplot(as.numeric(GM12878['IGKC',]), main = 'IGKC', ylab = 'Expression Level (CPM)', border = 'dodgerblue4', yaxt = 'n', ylim=c(0,35e4))
axis(2, at = seq(0, 300000, 300000/3), labels = formatC(seq(0, 300000, 300000/3), digits = 0, format = 'e'))
mtext('Cell Index', side = 1, cex = 0.7)
box()
barplot(as.numeric(GM12878['LTB',]), main = 'LTB', ylab = 'Expression Level (CPM)', border = 'dodgerblue4', yaxt = 'n', ylim = c(0,13000))
axis(2, at = seq(0,12000,12000/3), labels = formatC(seq(0, 12000,12000/3), digits = 0, format = 'e'))
mtext('Cell Index', side = 1, cex = 0.7)
box()
barplot(as.numeric(GM12878['RPL17',]), main = 'RPL17', ylab = 'Expression Level (CPM)', yaxt = 'n')
axis(2, at = seq(0, 1400, 1400/3), labels = formatC(seq(0, 1400, 1400/3), digits = 0, format = 'e'))
mtext('Cell Index', side = 1, cex = 0.7)
box()
barplot(as.numeric(GM12878['CCL3',]), main = 'CCL3', ylab = 'Expression Level (CPM)', yaxt = 'n',  border = 'dodgerblue4', ylim = c(0,27000))
axis(2, at = seq(0, 25000, 5000), labels = formatC(seq(0, 25000, 5000), digits = 0, format = 'e'))
mtext('Cell Index', side = 1, cex = 0.7)
box()
barplot(as.numeric(GM12878['FTL',]), main = 'FTL', ylab = 'Expression Level (CPM)', yaxt = 'n', border = 'dodgerblue4', ylim = c(0,64500))
axis(2, at = seq(0, 60000, 20000), labels = formatC(seq(0, 60000, 20000), digits = 0, format = 'e'))
mtext('Cell Index', side = 1, cex = 0.7)
box()
dev.off()
