library(Matrix)
library(RSpectra)
library(uwot)
library(igraph)
library(Seurat)
library(fgsea)

#### load data and packages
setwd("/data/victorxun/scTenifoldtime/GSE137524")
source("new_Normalization.R")
source("pcNet.R")
source("regression_fun.R")
source("tensorDecomposition.R")
source("UMAP_order.R")
source("scTenifoldTime.R")

dta <- readRDS("data/dta_raw.rds")
colnames(dta) <- gsub(pattern = "[.]", replacement = "_", x = colnames(dta))
dta_sudotime1 <- read.csv("data/Ordered_SCC6_repeat1.csv")
rownames(dta_sudotime1) <- gsub(pattern = "[.]", replacement = "_", x = dta_sudotime1$X)
dta_sudotime1 <- dta_sudotime1[, -1]
dta_sudotime1 <- dta_sudotime1[order(dta_sudotime1$O2), ]
dta <- dta[, rownames(dta_sudotime1)]
# dta <- scTenifoldNet::scQC(as.matrix(dta))
dta <- scTenifoldNet::scQC(as.matrix(dta), minPCT = 0.25)
dta <- new_Normalization(dta)
n_cell <- ncol(dta)
n_gene <- nrow(dta)

#### Now consider scTenifoldTime
dta_list <- list()
time_vec <- rep(NA, 9)
for (i in 1:8) {
  dta_list[[i]] <- as.matrix(dta[, (500 * (i - 1) + 1):(500 * (i + 1))])
  time_vec[i] <- mean(dta_sudotime1[colnames(dta_list[[i]]), 3])
}
dta_list[[9]] <- as.matrix(dta[, 4001:n_cell])
time_vec[9] <- mean(dta_sudotime1[colnames(dta_list[[9]]), 3])

set.seed(1)
res <- scTenifoldTime(dta_list = dta_list, time_vec = time_vec, nComp = 5, q = 0,
                                  K = 5, maxIter = 10000, maxError = 1e-5, thres = 0.05, nDecimal = 2,
                                  ma_nDim = 3)
saveRDS(res, "results_7000/pcnet_time1.rds")