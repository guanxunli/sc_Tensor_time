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
  network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
}

## community detection
# Getting beta mat
n_net <- length(network_list)
X_mat <- cbind(1, 1:n_net)
tmp <- as.numeric(network_list[[1]])
Y <- tmp
for (i in 2:n_net){
  tmp <- as.numeric(network_list[[i]])
  Y <- rbind(Y, tmp)
}
## Do regress
beta_coef <- solve(crossprod(X_mat), crossprod(X_mat, Y))
beta_mat <- matrix(beta_coef[2, ], nrow = nGenes)

res_APC <- APC_fun(beta_mat)
res_SCORE <- SCORE_fun(beta_adj, K = length(unique(res_APC)))

res <- list()
res$network <- network_list
res$beta_mat <- beta_mat
res$APC <- res_APC
res$SCORE <- res_SCORE
saveRDS(res, "res_cor.rds")

#### check results.
gene_list <- rownames(dta)
res <- readRDS("results/correlation method/res_cor.rds")
res_AP <- res$APC
res_SCORE <- res$SCORE
res_SLIM_AP <- res$SLIM_AP
names(res_AP) <- gene_list
names(res_SCORE) <- gene_list
names(res_SLIM_AP) <- gene_list

AP_list <- list()
i <- 1
tmp_community <- unique(res_AP)
for (iter in 1:length(unique(res_AP))) {
  tmp <- names(res_AP)[which(res_AP == tmp_community[iter])]
  if (length(tmp) > 5){
    AP_list[[i]] <- tmp
    i <- i + 1
  }
}

SCORE_list <- list()
i <- 1
tmp_community <- unique(res_SCORE)
for (iter in 1:length(unique(res_SCORE))) {
  tmp <- names(res_SCORE)[which(res_SCORE == tmp_community[iter])]
  if (length(tmp) > 5){
    SCORE_list[[i]] <- tmp
    i <- i + 1
  }
}

SLIM_AP_list <- list()
i <- 1
tmp_community <- unique(res_SLIM_AP)
for (iter in 1:length(unique(res_SLIM_AP))) {
  tmp <- names(res_SLIM_AP)[which(res_SLIM_AP == tmp_community[iter])]
  if (length(tmp) > 5){
    SLIM_AP_list[[i]] <- tmp
    i <- i + 1
  }
}

res <- list()
res$AP <- AP_list
res$SCORE <- SCORE_list
res$SLIM_AP <- SLIM_AP_list
saveRDS(res, "results/correlation method/res_community.rds")


