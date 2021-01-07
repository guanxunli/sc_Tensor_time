library(Matrix)
source("utility.R")
nCells <- c(10, 50, 100, 500, 1000, 2000, 3000)
## Z-normalization
metricOut <- list()
metricOut$Acc <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Rec <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Fmeasure <- matrix(NA, nrow = 10, ncol = length(nCells))
for (n in 1:length(nCells)){
  ncell <- nCells[n]
  for (i in 1:10){
    net_tmp <- readMM(paste0("networks/Z-normalization/Z-normalization_run", i, "_", ncell, "cells"))
    metricOut$Acc[i, n] <- Acurracy_fun(net_tmp, Q = 0)
    metricOut$Rec[i, n] <- Recall_fun(net_tmp, Q = 0)
    metricOut$Fmeasure[i, n] <- Fmeasure_fun(net_tmp, Q = 0)
  }
}

PCR_Znormalization <- metricOutput_fun(metricOut$Acc, metricOut$Rec, metricOut$Fmeasure, nCells = nCells)
PCR_Znormalization <- round(PCR_Znormalization, 3)

## log-normalization
metricOut <- list()
metricOut$Acc <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Rec <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Fmeasure <- matrix(NA, nrow = 10, ncol = length(nCells))
for (n in 1:length(nCells)){
  ncell <- nCells[n]
  for (i in 1:10){
    net_tmp <- readMM(paste0("networks/log-normalization/log-normalization_run", i, "_", ncell, "cells"))
    metricOut$Acc[i, n] <- Acurracy_fun(net_tmp, Q = 0)
    metricOut$Rec[i, n] <- Recall_fun(net_tmp, Q = 0)
    metricOut$Fmeasure[i, n] <- Fmeasure_fun(net_tmp, Q = 0)
  }
}

PCR_lognormalization <- metricOutput_fun(metricOut$Acc, metricOut$Rec,  metricOut$Fmeasure, nCells = nCells)
PCR_lognormalization <- round(PCR_lognormalization, 3)

## original normalization
metricOut <- list()
metricOut$Acc <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Rec <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Fmeasure <- matrix(NA, nrow = 10, ncol = length(nCells))
for (n in 1:length(nCells)){
  ncell <- nCells[n]
  for (i in 1:10){
    net_tmp <- readMM(paste0("networks/ori-normalization/ori-normalization_run", i, "_", ncell, "cells"))
    metricOut$Acc[i, n] <- Acurracy_fun(net_tmp, Q = 0)
    metricOut$Rec[i, n] <- Recall_fun(net_tmp, Q = 0)
    metricOut$Fmeasure[i, n] <- Fmeasure_fun(net_tmp, Q = 0)
  }
}

PCR_orinormalization <- metricOutput_fun(metricOut$Acc, metricOut$Rec, metricOut$Fmeasure, nCells = nCells)
PCR_orinormalization <- round(PCR_orinormalization, 3)

## new normalization
metricOut <- list()
metricOut$Acc <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Rec <- matrix(NA, nrow = 10, ncol = length(nCells))
metricOut$Fmeasure <- matrix(NA, nrow = 10, ncol = length(nCells))
for (n in 1:length(nCells)){
  ncell <- nCells[n]
  for (i in 1:10){
    net_tmp <- readMM(paste0("networks/new-normalization/new-normalization_run", i, "_", ncell, "cells"))
    metricOut$Acc[i, n] <- Acurracy_fun(net_tmp, Q = 0)
    metricOut$Rec[i, n] <- Recall_fun(net_tmp, Q = 0)
    metricOut$Fmeasure[i, n] <- Fmeasure_fun(net_tmp, Q = 0)
  }
}

PCR_newnormalization <- metricOutput_fun(metricOut$Acc, metricOut$Rec, metricOut$Fmeasure, nCells = nCells)
PCR_newnormalization <- round(PCR_newnormalization, 3)

## plot results
library(RColorBrewer)
pColors <- RColorBrewer::brewer.pal(4,'Dark2')

png('figures/ACC.png', width = 900, height = 900, res = 300, pointsize = 10)
par(mar=c(3.5,3.5,1,1), mgp=c(2.4,0.5,0), las = 1)
plot(PCR_Znormalization$nCells, PCR_Znormalization$acc, ylim=c(0.5,0.9), type='b', col = pColors[1],
     ylab = '', xlab = '', pch = 15)
mtext('Number of Cells', side = 1, line = 1.6)
mtext('Accuracy', side = 2, line = 2.5, las = 3)
arrows(x0 = PCR_Znormalization$nCells, x1 = PCR_Znormalization$nCells, y0 = PCR_Znormalization$accLB, y1 = PCR_Znormalization$accUB,
       length = 0.03, code = 3, angle = 90, col = pColors[1])
points(PCR_lognormalization$nCells, PCR_lognormalization$acc, col = pColors[2], type='b', pch = 16)
arrows(x0 = PCR_lognormalization$nCells, x1 = PCR_lognormalization$nCells, y0 = PCR_lognormalization$accLB, y1 = PCR_lognormalization$accUB,
       length = 0.03, code = 3, angle = 90, col = pColors[2])
points(PCR_orinormalization$nCells, PCR_orinormalization$acc, col = pColors[3], type='b', pch = 17)
arrows(x0 = PCR_orinormalization$nCells, x1 = PCR_orinormalization$nCells, y0 = PCR_orinormalization$accLB, y1 = PCR_orinormalization$accUB,
       length = 0.03, code = 3, angle = 90, col = pColors[3])
points(PCR_newnormalization$nCells, PCR_newnormalization$acc, col = pColors[4], type='b', pch = 18)
arrows(x0 = PCR_newnormalization$nCells, x1 = PCR_newnormalization$nCells, y0 = PCR_newnormalization$accLB, y1 = PCR_newnormalization$accUB,
       length = 0.03, code = 3, angle = 90, col = pColors[4])
legend('bottomright', legend = c('PCR_Z','PCR_log', 'PCR_ori', 'PCR_new'), bty = 'n',
       cex = 0.7, ncol = 2, pch = 15:18, col = pColors)
dev.off()

png('figures/REC.png', width = 900, height = 900, res = 300, pointsize = 10)
par(mar=c(3.5,3.5,1,1), mgp=c(2.4,0.5,0), las = 1)
plot(PCR_Znormalization$nCells, PCR_Znormalization$recall, ylim=c(0.15,1), type='b', col = pColors[1],
     ylab = '', xlab = '', pch = 15)
mtext('Number of Cells', side = 1, line = 1.6)
mtext('Recall', side = 2, line = 2.5, las = 3)
arrows(x0 = PCR_Znormalization$nCells, x1 = PCR_Znormalization$nCells, y0 = PCR_Znormalization$recallLB, y1 = PCR_Znormalization$recallUB,
       length = 0.03, code = 3, angle = 90, col = pColors[1])
points(PCR_lognormalization$nCells, PCR_lognormalization$recall, col = pColors[2], type='b', pch = 16)
arrows(x0 = PCR_lognormalization$nCells, x1 = PCR_lognormalization$nCells, y0 = PCR_lognormalization$recallLB, y1 = PCR_lognormalization$recallUB,
       length = 0.03, code = 3, angle = 90, col = pColors[2])
points(PCR_orinormalization$nCells, PCR_orinormalization$recall, col = pColors[3], type='b', pch = 17)
arrows(x0 = PCR_orinormalization$nCells, x1 = PCR_orinormalization$nCells, y0 = PCR_orinormalization$recallLB, y1 = PCR_orinormalization$recallUB,
       length = 0.03, code = 3, angle = 90, col = pColors[3])
points(PCR_newnormalization$nCells, PCR_newnormalization$recall, col = pColors[4], type='b', pch = 18)
arrows(x0 = PCR_newnormalization$nCells, x1 = PCR_newnormalization$nCells, y0 = PCR_newnormalization$recallLB, y1 = PCR_newnormalization$recallUB,
       length = 0.03, code = 3, angle = 90, col = pColors[4])
legend('bottomright', legend = c('PCR_Z','PCR_log', 'PCR_ori', 'PCR_new'), bty = 'n',
       cex = 0.7, ncol = 2, pch = 15:18, col = pColors)
dev.off()

png('figures/Fmeasure.png', width = 900, height = 900, res = 300, pointsize = 10)
par(mar=c(3.5,3.5,1,1), mgp=c(2.4,0.5,0), las = 1)
plot(PCR_Znormalization$nCells, PCR_Znormalization$Fmeasure, ylim=c(0.15,1), type='b', col = pColors[1],
     ylab = '', xlab = '', pch = 15)
mtext('Number of Cells', side = 1, line = 1.6)
mtext('Fmeasure', side = 2, line = 2.5, las = 3)
arrows(x0 = PCR_Znormalization$nCells, x1 = PCR_Znormalization$nCells, y0 = PCR_Znormalization$FmeasureLB, y1 = PCR_Znormalization$FmeasureUB,
       length = 0.03, code = 3, angle = 90, col = pColors[1])
points(PCR_lognormalization$nCells, PCR_lognormalization$Fmeasure, col = pColors[2], type='b', pch = 16)
arrows(x0 = PCR_lognormalization$nCells, x1 = PCR_lognormalization$nCells, y0 = PCR_lognormalization$FmeasureLB, y1 = PCR_lognormalization$FmeasureUB,
       length = 0.03, code = 3, angle = 90, col = pColors[2])
points(PCR_orinormalization$nCells, PCR_orinormalization$Fmeasure, col = pColors[3], type='b', pch = 17)
arrows(x0 = PCR_orinormalization$nCells, x1 = PCR_orinormalization$nCells, y0 = PCR_orinormalization$FmeasureLB, y1 = PCR_orinormalization$FmeasureUB,
       length = 0.03, code = 3, angle = 90, col = pColors[3])
points(PCR_newnormalization$nCells, PCR_newnormalization$Fmeasure, col = pColors[4], type='b', pch = 18)
arrows(x0 = PCR_newnormalization$nCells, x1 = PCR_newnormalization$nCells, y0 = PCR_newnormalization$FmeasureLB, y1 = PCR_newnormalization$FmeasureUB,
       length = 0.03, code = 3, angle = 90, col = pColors[4])
legend('bottomright', legend = c('PCR_Z','PCR_log', 'PCR_ori', 'PCR_new'), bty = 'n',
       cex = 0.7, ncol = 2, pch = 15:18, col = pColors)
dev.off()