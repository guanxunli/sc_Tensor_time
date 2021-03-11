library(Matrix)
library(RSpectra)

#### load data and packages
setwd("/data/victorxun/scTenifoldtime/zebrafish_hindbrain")
source("new_Normalization.R")
source("pcNet.R")
source("regression_fun.R")
source("tensorDecomposition.R")
source("UMAP_order.R")
source("scTenifoldTime.R")

######## For group 2
dta <- readRDS("data/dta_small.rds")
sudotime <- readRDS("data/sudotime_small.rds")
dta <- scTenifoldNet::scQC(as.matrix(dta)) 
dta <- new_Normalization(dta)
n_cell <- ncol(dta)
n_gene <- nrow(dta)
res <- list()

#### Now consider scTenifoldTime
dta_list <- list()
time_vec <- rep(NA, 5)
for (i in 1:4) {
  dta_list[[i]] <- as.matrix(dta[, (250 * (i - 1) + 1):(250 * (i + 1))])
  time_vec[i] <- mean(sudotime[colnames(dta_list[[i]]), 1])
}
dta_list[[5]] <- as.matrix(dta[, 1000:n_cell])
time_vec[5] <- mean(sudotime[colnames(dta_list[[5]]), 1])
time_vec <- time_vec / max(time_vec)

#### correlation method
res_small <- scTenifoldTime_beta(dta_list, method = "pcnet", time_vec, nComp = 5, q = 0,
                                 K = 5, maxIter = 10000, maxError = 1e-5)
res$res_small <- res_small
rm(res_small)

######### For all cells
dta <- readRDS("data/dta_large.rds")
sudotime <- readRDS("data/sudotime_large.rds")
dta <- scTenifoldNet::scQC(as.matrix(dta)) 
dta <- new_Normalization(dta)
n_cell <- ncol(dta)
n_gene <- nrow(dta)

#### Now consider scTenifoldTime
dta_list <- list()
time_vec <- rep(NA, 10)
for (i in 1:9) {
  dta_list[[i]] <- as.matrix(dta[, (500 * (i - 1) + 1):(500 * (i + 1))])
  time_vec[i] <- mean(sudotime[colnames(dta_list[[i]]), 1])
}
dta_list[[10]] <- as.matrix(dta[, 4500:n_cell])
time_vec[10] <- mean(sudotime[colnames(dta_list[[10]]), 1])
time_vec <- time_vec / max(time_vec)

#### correlation method
res_large<- scTenifoldTime_beta(dta_list, method = "pcnet", time_vec, nComp = 5, q = 0,
                                K = 5, maxIter = 10000, maxError = 1e-5)
res$res_large <- res_large
rm(res_large)

saveRDS(res, "results/beta_pcnet.rds")