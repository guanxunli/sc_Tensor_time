library(Matrix)
library(RSpectra)

#### load data and packages
setwd("/data/victorxun/scTenifoldtime/GSE147807")
source("new_Normalization.R")
source("pcNet.R")
source("regression_fun.R")
source("tensorDecomposition.R")
source("UMAP_order.R")
source("scTenifoldTime.R")

dta <- readRDS("data/dta_raw.rds")
dta_sudotime <- readRDS("data/dta_sudotime.rds")
dta <- scTenifoldNet::scQC(as.matrix(dta),  maxMTratio = 0.5, minPCT = 0.25) # Is it OK for maxMTratio = 0.5
dta <- new_Normalization(dta)
n_cell <- ncol(dta)
n_gene <- nrow(dta)
res <- list()

#### Now consider scTenifoldTime
dta_list <- list()
time_vec <- rep(NA, 5)
for (i in 1:4) {
  dta_list[[i]] <- as.matrix(dta[, (250 * (i - 1) + 1):(250 * (i + 1))])
  time_vec[i] <- mean(dta_sudotime[colnames(dta_list[[i]]), 1])
}
dta_list[[5]] <- as.matrix(dta[, 1000:n_cell])
time_vec[5] <- mean(dta_sudotime[colnames(dta_list[[5]]), 1])

#### correlation method
res_cor <- scTenifoldTime_beta(dta_list, method = "pcnet", time_vec, nComp = 5, q = 0,
                               K = 10, maxIter = 10000, maxError = 1e-5)

saveRDS(res, "results/beta_pcnet.rds")