library(Matrix)
source("utility.R")

## load data
dta <- readRDS("data/dta.rds")
dta_sudotime <- read.csv("data/monocle_rank.csv")
dta_sudotime <- dta_sudotime[order(dta_sudotime$CCAT_score), ]
dta <- dta[, dta_sudotime$Barcode]
dta <- scTenifoldNet::scQC(as.matrix(dta))
dta <- new_Normalization(dta)

dta_list <- list()
for (i in 1:6){
  dta_list[[i]] <- as.matrix(dta[, (400 * (i - 1) + 1):(400 * (i + 1))])
}
dta_list[[7]] <- as.matrix(dta[, 2401:3120])

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
beta_mat <- matrix(beta_coef[2, ], nrow =  nrow = nrow(network_list[[1]]))
beta_mat_abs <- round(abs(beta_mat), 2)
beta_adj <- beta_mat_abs
beta_adj[which(beta_adj > 0)] = 1

res_APC <- APC_fun(beta_mat)
res_SCORE <- SCORE_fun(beta_adj, K = length(unique(res_APC)))
res_SLIM <- MASE_fun(network_list, di = 50, K = 50, gamma = 0.25, m = 8)

res <- list()
res$network <- network_list
res$APC <- res_APC
res$SCORE <- res_SCORE
res$SLIM_AP <- res_SLIM
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

#### Check community results
res <- readRDS("results/correlation method/res_community.rds")
library(enrichR)
library(igraph)

score_function <- function(res_list, fdrThreshold = 0.05){
  n <- length(res_list)
  n_community <- 0
  for (i in 1:length(res_list)){
    gList <- res_list[[i]]
    gList <- gList[!grepl('^mt-|^Rpl|^Rps',gList, ignore.case = TRUE)]
    if (length(gList) == 0){
      n_community <- n_community
    } else{
      E <- enrichr(gList, c('KEGG_2019_Mouse','Reactome_2016','WikiPathways_2019_Mouse'))
      E <- do.call(rbind.data.frame, E)
      E <- E[E$Adjusted.P.value < fdrThreshold,]
      if(isTRUE(nrow(E) > 0)){
        n_community <- n_community  + 1
      }
    }
  }
  return(n_community / n)
}

res_score <- lapply(res, score_function)
res_score
