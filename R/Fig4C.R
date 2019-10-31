library(Matrix)
library(phateR)
library(Seurat)
library(plot3D)
library(rgl)

# Reading imputed values
igeneMatrix <- t(read.csv('../Data/nfkb8genes.csv'))
geneList <- rownames(igeneMatrix)

# Associated Values
nfkbAssociations <- read.csv('../Data/nfkbAssociations.csv', row.names = 1)

png('../Results/figures/FIG4C.png', width = 2000, height = 2000, res = 300, pointsize = 20)
par(mfrow=c(8,8), mar=c(.5,.5,.5,.5))
for(i in 1:8){
  for(j in 1:8){
    if(i == j){
        hist((igeneMatrix[i,]), main = '',xaxt='n', yaxt='n')
      if(i == 1 & j == 1){
        mtext(text = geneList[j], line = -0.1, side = 3, cex = 0.4)
        mtext(text = geneList[i], line = -0.1, side = 2, cex = 0.4)
      }
    } else {
      if(i == 1){
        data <- cbind(igeneMatrix[i,], igeneMatrix[j,])
        data <- densCols(data)
        if(!is.na(nfkbAssociations[i,j])){
        plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = data)
        abline(lm(igeneMatrix[j,]~igeneMatrix[i,]), col='red')
        corVal <- cor(igeneMatrix[i,], igeneMatrix[j,], method = 'spearman')
        legend('topleft', legend = parse(text = paste0('rho ==', round(corVal,2))), bty = 'n', cex = 0.4)
        if(nfkbAssociations[i,j] > 0 & corVal> 0){
            box(col='forestgreen', lwd=3)  
          }
          if(nfkbAssociations[i,j] <  0 & corVal <  0){
            box(col='forestgreen', lwd=3)  
          }
          if(nfkbAssociations[i,j] <  0 & corVal >  0){
            box(col='red', lwd=3)  
          }
          if(nfkbAssociations[i,j] >  0 & corVal <  0){
            box(col='red', lwd=3)  
          }
        } else {
          plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = 'gray90')
          box(col='gray90', lwd=3)
        }
        mtext(text = geneList[j], line = -0.1, side = 3, cex = 0.4)
        
      }else{
        if(j == 1){
          data <- cbind(igeneMatrix[i,], igeneMatrix[j,])
          data <- densCols(data)
          if(!is.na(nfkbAssociations[i,j])){
            plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = data)
            abline(lm(igeneMatrix[j,]~igeneMatrix[i,]), col='red')
            corVal <- cor(igeneMatrix[i,], igeneMatrix[j,], method = 'spearman')
            legend('topleft', legend = parse(text = paste0('rho ==', round(corVal,2))), bty = 'n', cex = 0.4)
            if(nfkbAssociations[i,j] > 0 & corVal> 0){
              box(col='forestgreen', lwd=3)  
            }
            if(nfkbAssociations[i,j] <  0 & corVal <  0){
              box(col='forestgreen', lwd=3)  
            }
            if(nfkbAssociations[i,j] <  0 & corVal >  0){
              box(col='red', lwd=3)  
            }
            if(nfkbAssociations[i,j] >  0 & corVal <  0){
              box(col='red', lwd=3)  
            }
          } else {
            plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = 'gray90')
            box(col='gray90', lwd=3)
          }
          mtext(text = geneList[i], line = -0.1, side = 2, cex = 0.4)
        }else {
          data <- cbind(igeneMatrix[i,], igeneMatrix[j,])
          data <- densCols(data)
          if(!is.na(nfkbAssociations[i,j])){
            plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = data)
            abline(lm(igeneMatrix[j,]~igeneMatrix[i,]), col='red')
            corVal <- cor(igeneMatrix[i,], igeneMatrix[j,], method = 'spearman')
            legend('topleft', legend = parse(text = paste0('rho ==', round(corVal,2))), bty = 'n', cex = 0.4)
            if(nfkbAssociations[i,j] > 0 & corVal> 0){
              box(col='forestgreen', lwd=3)  
            }
            if(nfkbAssociations[i,j] <  0 & corVal <  0){
              box(col='forestgreen', lwd=3)  
            }
            if(nfkbAssociations[i,j] <  0 & corVal >  0){
              box(col='red', lwd=3)  
            }
            if(nfkbAssociations[i,j] >  0 & corVal <  0){
              box(col='red', lwd=3)  
            }
          } else {
            plot(igeneMatrix[i,], igeneMatrix[j,], cex = 0.1, pch = 16, xaxt='n', yaxt='n', col = 'gray90')
            box(col='gray90', lwd=3)
          }
        }
      }
    }
  }
}
dev.off()

