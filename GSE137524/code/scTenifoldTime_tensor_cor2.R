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
dta_sudotime2 <- read.csv("data/Ordered_SCC6_repeat2.csv")
dta_sudotime2 <- dta_sudotime2[order(dta_sudotime2$O2), ]
rownames(dta_sudotime2) <- gsub(pattern = "[.]", replacement = "_", x = dta_sudotime2$X)
dta_sudotime2 <- dta_sudotime2[, -1]
dta <- dta[, rownames(dta_sudotime2)]
# dta <- scTenifoldNet::scQC(as.matrix(dta))
dta <- scTenifoldNet::scQC(as.matrix(dta), minPCT = 0.25)
dta <- new_Normalization(dta)
n_cell <- ncol(dta)
n_gene <- nrow(dta)

#### Now consider our own method
dta_list <- list()
time_vec <- rep(NA, 7)
for (i in 1:6) {
  dta_list[[i]] <- as.matrix(dta[, (400 * (i - 1) + 1):(400 * (i + 1))])
  time_vec[i] <- mean(dta_sudotime2[colnames(dta_list[[i]]), 3])
}

dta_list[[7]] <- as.matrix(dta[, 2401:n_cell])
time_vec[7] <- mean(dta_sudotime2[colnames(dta_list[[7]]), 3])
time_vec <- time_vec / max(time_vec)


# dta_list <- list()
# time_vec <- rep(NA, 7)
# for (i in 1:6) {
#   dta_list[[i]] <- as.matrix(dta[, (500 * (i - 1) + 1):(500 * (i + 1))])
#   time_vec[i] <- mean(dta_sudotime2[colnames(dta_list[[i]]), 3])
# }
# dta_list[[7]] <- as.matrix(dta[, 2501:n_cell])
# time_vec[7] <- mean(dta_sudotime2[colnames(dta_list[[7]]), 3])
# 
set.seed(1)
res <- scTenifoldTime_tensor(dta_list = dta_list, time_vec = time_vec, method = "cor", nComp = 5, q = 0,
                      K = 10, maxIter = 10000, maxError = 1e-5, thres = 0.05, nDecimal = 2,
                      ma_nDim = 30)

saveRDS(res, "results_7000/tensor_cor2.rds")
# saveRDS(res, "results_10000/pcnet_time2.rds")