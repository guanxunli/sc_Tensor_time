library(Seurat)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/filterGenes.R')
setwd("Z:/Cailab/pseudotime_dataset/GSE137524_HNSCC_cells")
SCC6 <- read.csv('GSE137524_exprsSCC6Matrix.csv')
rownames(SCC6) <- make.unique(SCC6$X)
SCC6 <- SCC6[,-c(1)]
SCC6 <-scQC(SCC6)
SCC6 <- filterGenes(SCC6)

cellnumber <- colnames(SCC6)
SCC6_CTX_R1ID <- which(grepl('SCC6_CTX_R1',cellnumber))
SCC6_PBS_R1ID <- which(grepl('SCC6_PBS_R1',cellnumber))

CTX <- SCC6[,c(SCC6_CTX_R1ID)]
PBS <- SCC6[,c(SCC6_PBS_R1ID)]
CTX <- CreateSeuratObject(count=CTX,project = "CTX")
PBS <- CreateSeuratObject(count=PBS,project = "PBS")
A <- merge(CTX,PBS)
table(Idents(A))
A <- FindVariableFeatures(A)
A <- ScaleData(A)
A <- RunPCA(A)
A <- RunUMAP(A, dims = 1:20)
UMAPPlot(A)

B <- A@assays$RNA@counts
B <- B[VariableFeatures(A),]
suppressMessages(library(monocle))
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/getPseudotime.R')
library(monocle)
O2 <- computePseudoTime(data.matrix(B),nDim = 100,simplified = FALSE)

plot(attr(O2, 'DDRTree'))
library(phateR)
O1 <- phate(as.matrix(t(B)))
plot(O1$embedding)
P1 <- data.frame(O1$embedding, R=rank(O2))
P2 <- data.frame(attr(O2, 'DDRTree'), R = rank(O2))
P3 <- data.frame(A@reductions$umap@cell.embeddings, R=rank(O2))
P4 <- data.frame(A@reductions$pca@cell.embeddings, R=rank(O2))

pA <- ggplot(P2, aes(DDRTree1, DDRTree2, color=R)) + geom_point() + theme_bw() + labs(title = 'MONOCLE3')
pB <- ggplot(P1, aes(PHATE1, PHATE2, color=R)) + geom_point() + theme_bw() + labs(title = 'PHATE')
pC <- ggplot(P3, aes(UMAP_1, UMAP_2, color=R)) + geom_point() + theme_bw() + labs(title = 'UMAP')
pD <- ggplot(P4, aes(PC_1,PC_2,color=R)) + geom_point() + theme_bw() + labs(title = 'PCA')

setwd("Z:/Cailab/pseudotime_dataset/GSE137524_HNSCC_cells/OrganizeCells")
png('Order_SCC6_repeat1.png',width = 3500,height = 3000,res=300)
pA+pB+pC+pD
dev.off()
png('Umap_SCC6_repeat1.png',width = 1500,height = 1000,res=300)
UMAPPlot(A)
dev.off()

SCC6_repeat1 <- as.data.frame(cbind(attr(O2,"DDRTree"),O2))
CTXid <- which(grepl('SCC6_CTX_R1',rownames(SCC6_repeat1)))
PBSid <- which(grepl('SCC6_PBS_R1',rownames(SCC6_repeat1)))
SCC6_repeat1$treatment <- 'NA'
SCC6_repeat1$treatment[CTXid] <- 'CTX'
SCC6_repeat1$treatment[PBSid] <- 'PBS'
write.csv(SCC6_repeat1,file = 'Ordered_SCC6_repeat1.csv')

library(ggpubr)
library(ggplot2)
library(ggstatsplot)
png('Boxplot_SCC6_repeat1.png',width = 1600,height = 1200,res=300)
ggboxplot(data=SCC6_repeat1, x='treatment', y='O2', fill='treatment',
          add = "jitter",
          shape = "treatment",
          xlab = "Treatment",
          ylab = 'Pseudotime') + stat_compare_means(method = "t.test",
                                                    label.x = 1.35, 
                                                    label.y = 47)
dev.off()


rm(list=setdiff(ls(), "SCC6"))
