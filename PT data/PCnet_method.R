library(Matrix)
source("utility.R")

## load data
dta <- readRDS("data/dta.rds")
dta_sudotime <- read.csv("data/monocle_rank.csv")
dta_sudotime <- dta_sudotime[order(dta_sudotime$CCAT_score), ]
dta <- dta[, dta_sudotime$Barcode]
dta <- scTenifoldNet::scQC(dta)
dta <- new_Normalization(dta)

dta_list <- list()
for (i in 1:6){
  dta_list[[i]] <- dta[, (400 * (i - 1) + 1):(400 * (i + 1))]
}
dta_list[[7]] <- dta[, 2401:3120]

## make networks
network_list <- list()
for (i in 1:length(dta_list)){
  network_list[[i]] <- scTeni_makeNet(dta_list[[i]], norm_method = "new", nc_nNet = 10, nc_nCells = 500, nc_nComp = 5, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                                                  nc_q = 0.5, td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5)
}

# Getting beta mat
n_net <- length(network_list)
X_mat <- cbind(1, 1:n_net)
tmp <- as.numeric(network_list[[1]])
Y <- tmp
for (i in 2:10){
  tmp <- as.numeric(network_list[[i]])
  Y <- rbind(Y, tmp)
}
## Do regress
beta_coef <- solve(crossprod(X_mat), crossprod(X_mat, Y))
beta_mat <- matrix(beta_coef[2, ], nrow = nrow(network_list[[1]]))
res_AP <- APC_fun(beta_mat)

res_SLIM <- MASE_fun(network_list, di = 10, K = 10, gamma = 0.25, m = 8)
res <- list()
res$network <- network_list
res$beta_mat <- beta_mat
res$AP <- res_AP
res$SLIM <- res_SLIM
saveRDS(res, "res_PCnet400.rds")




