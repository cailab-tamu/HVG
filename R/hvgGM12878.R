library(Matrix)
cMatrix <- readMM('hvgGM12878/matrix.mtx')
colnames(cMatrix) <- readLines('hvgGM12878/barcodes.tsv')
rownames(cMatrix) <- read.csv('hvgGM12878/genes.tsv', sep = '\t', stringsAsFactors = FALSE, header = FALSE)[,2]
cMatrix <- t(t(cMatrix)/colSums(cMatrix)) * 1e6
source('https://raw.githubusercontent.com/dosorio/utilities/master/LMA/knnSearch.R')

#QC
mtRate <- colSums(cMatrix[grepl('MT-',rownames(cMatrix)),])/colSums(cMatrix)
cMatrix <- cMatrix[,mtRate < 0.1]
cMatrix <- t(t(cMatrix)/colSums(cMatrix)) * 1e6

# immuneGenes <- read.csv('GO_term_summary_20191022_094850.txt', sep = '\t', row.names = NULL, stringsAsFactors = FALSE)
# immuneGenes <- toupper(unique(immuneGenes$MGI.Gene.Marker.ID))
# cMatrix <- cMatrix[!rownames(cMatrix) %in% immuneGenes,]

source('https://raw.githubusercontent.com/dosorio/utilities/master/LMA/knnSearch.R')

k <- knnSearch(cMatrix, 1000)
rownames(k) <- colnames(k) <- NULL
k <- reshape2::melt(t(k))[,-2]
k <- igraph::graph_from_data_frame(k, directed = TRUE)
k <- k[]


set.seed(1)
nBoot <- 1000
sCells <- sample(seq_len(ncol(cMatrix)),nBoot, replace = TRUE)

HVG <- pbapply::pbsapply(sCells, function(X){
  temp <- cMatrix[,which(as.logical(k[X,]))]
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
  
  
  pVal <- pchisq(log(cv2/pFit)*1000,999, lower.tail = FALSE)
  pAdj <- p.adjust(pVal, method = 'fdr')
  #plot(xplot, col = ifelse(rownames(temp) %in% names(pAdj[pAdj < 0.05]),'red','black'))
  names(pAdj[pAdj < 0.05])
})

HVG <- sort(table(unlist(HVG)), decreasing = TRUE)
HVG2 <- names(HVG[(HVG/nBoot) >= 0.9])
write.csv(cbind(HVG/nBoot), 'LCL_HVG.csv')
HVG <- names(HVG)#[HVG >= quantile(HVG, 0.9)])

source('https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/hsa_GO_SYMBOL.R')
GO <- hsa_GO_SYMBOL(HVG)
GO2 <- hsa_GO_SYMBOL(HVG2)

GO <- GO[order(GO$p.adjust,decreasing = FALSE),]
GO2 <- GO2[order(GO2$p.adjust,decreasing = FALSE),]

write.csv(GO, file = 'LCL_GO_UNION.csv')
write.csv(GO2, file = 'LCL_GO_INTERSECTION.csv')

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

library(ggplot2)
library(ggrepel)
library(ggpubr)
png('LCL.png', width = 8000, height = 3000, res = 300)
gCol <- ifelse(rownames(xplot) %in% HVG, 'red', 'black')
p1 <- ggplot(xplot, mapping = aes(x = mean_log, y = cv2_log)) + 
  geom_text_repel(data = xplot[HVG,], label=HVG)+
  geom_point(color=gCol) + 
  theme_classic() + 
  xlab(parse(text = 'log[10](Mean)')) +
  ylab(parse(text = 'log[10](CV^2)')) +
  ggtitle(label = 'LYMPHOBLASTOID CELL LINE', subtitle = '(1000 randomizations) Union')
gCol <- ifelse(rownames(xplot) %in% HVG2, 'red', 'black')
p2 <- ggplot(xplot, mapping = aes(x = mean_log, y = cv2_log)) + 
  geom_text_repel(data = xplot[HVG2,], label=HVG2)+
  geom_point(color=gCol) + 
  theme_classic() + 
  xlab(parse(text = 'log[10](Mean)')) +
  ylab(parse(text = 'log[10](CV^2)')) +
  ggtitle(label = 'LYMPHOBLASTOID CELL LINE',subtitle = '(1000 randomizations) Intersection')
ggarrange(p1,p2)
dev.off()