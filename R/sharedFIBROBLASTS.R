library(VennDiagram)
DF <- read.csv('../Results/hvgList/DF_HVG.csv', stringsAsFactors = FALSE, row.names = 1)
DF <- DF$GENE[DF$FREQUENCY == 1]

FIBROBLAST1 <- read.csv('../Results/hvgList/FIBROBLAST1_HVG.csv', stringsAsFactors = FALSE, row.names = 1)
FIBROBLAST1 <- FIBROBLAST1$GENE[FIBROBLAST1$FREQUENCY == 1]

FIBROBLAST2 <- read.csv('../Results/hvgList/FIBROBLAST2_HVG.csv', stringsAsFactors = FALSE, row.names = 1)
FIBROBLAST2 <- FIBROBLAST2$GENE[FIBROBLAST2$FREQUENCY == 1]

sharedGenes <- list(DF,FIBROBLAST1,FIBROBLAST2)
names(sharedGenes) <- c('Dermal\nFibroblasts', 'Lung Distal\nFibroblasts', 'Lung Proximal\nFibroblasts')

VennDiagram::venn.diagram(x = sharedGenes, filename = '../Results/figures/sharedFIBROBLASTS.png', height = 1500, width = 1500, resolution = 300,
                          fill=c('blue','red','green'), cat.fontfamily='sans',sep.dist = 0.1, rotation.degree = 30,lty = "blank",cat.pos = c(-45,35,135),
                          fontfamily= 'sans')

gList <- table(unlist(sharedGenes))
sList <- (names(gList[gList == 3]))

source('https://raw.githubusercontent.com/dosorio/utilities/master/enrichments/hsa_GO_SYMBOL.R')
sEnrichment <- hsa_GO_SYMBOL(sList)
sEnrichment <- sEnrichment[order(sEnrichment$p.adjust, decreasing = FALSE),]

write.csv(sEnrichment, '../Results/annotations/sharedFIBROBLASTS.csv')
