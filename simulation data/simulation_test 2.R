library(Matrix)
source("simulation_utility.R")
## load data
dta <- read.csv("SERGIO/diffOutput_S_1000.csv", header = FALSE)
rownames(dta) <- paste0("Gene", 1:nrow(dta))
colnames(dta) <- paste0("Cell", 1:ncol(dta))
nGenes <- nrow(dta)

dta_list <- list()
for (i in 1:10){
  dta_list[[i]] <- dta[, (1000 * (i - 1) + 1):(1000 * i)]
}

npc_index <- c(3, 5, 9, 13)
ncell_index <- c(500, 700, 800, 900)
nq_index <- c(0.05, 0.25, 0.5)
tensor_index <- c(0, 3, 20)
community_index <- c("AP", "SCORE", "LPA", "SLIM", "M_SCORE", "M_SLIM", "M_SCORE_AP", "M_SCORE_LPA", "M_SLIM_AP", "M_SLIM_LPA")
index <- expand.grid(npc_index = npc_index, ncell_index = ncell_index,
                     tensor_index = tensor_index, nq_index = nq_index, community_index = community_index)

## define function
network_fun <- function(dta_list, method_index, npc_index, ncell_index, tensor_index, nq_index, community_index){
  ## make networks
  network_list <- list()
  for (i in 1:length(dta_list)){
    if (method_index == "cor"){
      network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
    } else{
      network_list[[i]] <- scTeni_makeNet(dta_list[[i]], norm_method = "old", nc_nNet = 10, nc_nCells = ncell_index, nc_nComp = npc_index, 
                                          nc_q = nq_index, td_K = tensor_index)
    }
  }
  
  ## community detection
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
  beta_mat <- matrix(beta_coef[2, ], nrow = nGenes)
  beta_mat_abs <- abs(beta_mat)
  beta_adj <- round(beta_mat_abs, 2)
  beta_adj[beta_adj > 0] <- 1
  # Do community detection
  true_cluster <- c(rep(1, 29), rep(2, 21), rep(3, 30), rep(4, 20), 1, 3, 4)
  if (community_index == "AP"){
    res <- APC_fun(beta_mat)
  } else if(community_index == "LPA"){
    res <- LPA_fun(beta_mat, K = 5)
  } else if(community_index == "SCORE"){
    res <- SCORE_fun(beta_adj, K = 4)
  } else if (community_index == "SLIM"){
    res <- SLIM_fun(beta_mat_abs, K = 4)
  } else if (community_index == "M_SCORE"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SCORE", c_method = "kmeans")
  } else if (community_index == "M_SLIM"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SLIM", c_method = "kmeans")
  } else if (community_index == "M_SCORE_AP"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SCORE", c_method = "AP")
  } else if (community_index == "M_SCORE_LPA"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SCORE", c_method = "LPA")
  } else if (community_index == "M_SLIM_AP") {
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SLIM", c_method = "AP")
  } else if (community_index == "M_SLIM_LPA"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SLIM", c_method = "LPA")
  }
  # return results
  return(fossil::rand.index(res, true_cluster))
}

## Do simulation
library(foreach)
library(doParallel)
library(doRNG)

## cor method
randindex_cor <- rep(NA, length(community_index))
for (i in 1:length(community_index)){
  randindex_cor[i] <- network_fun(dta_list, method_index = "cor", community_index = community_index[i])
}
saveRDS(randindex_cor, "out_cor_withoutnorm1000.rds")

cl <- makeCluster(20) # 4 workers
registerDoParallel(cl)

set.seed(0104)
out_res <- foreach(iter = 1:nrow(index)) %dorng% {
  try(network_fun(dta_list, method_index = "PCnet", npc_index = index[iter, 1], ncell_index = index[iter, 2], 
                  tensor_index = index[iter, 3], nq_index = index[iter, 4], community_index = index[iter, 5]))
}
stopCluster(cl)
saveRDS(out_res, "out_res_withoutnorm_1000.rds")

########################################## With new normalization method ##########################################
dta <- new_Normalization(dta)
dta_list <- list()
for (i in 1:10){
  dta_list[[i]] <- dta[, (1000 * (i - 1) + 1):(1000 * i)]
}

## define function
network_fun <- function(dta_list, method_index, npc_index, ncell_index, tensor_index, nq_index, community_index){
  ## make networks
  network_list <- list()
  for (i in 1:length(dta_list)){
    if (method_index == "cor"){
      network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
    } else{
      network_list[[i]] <- scTeni_makeNet(dta_list[[i]], norm_method = "new", nc_nNet = 10, nc_nCells = ncell_index, nc_nComp = npc_index, 
                                          nc_q = nq_index, td_K = tensor_index)
    }
  }
  
  ## community detection
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
  beta_mat <- matrix(beta_coef[2, ], nrow = nGenes)
  beta_mat_abs <- abs(beta_mat)
  beta_adj <- round(beta_mat_abs, 2)
  beta_adj[beta_adj > 0] <- 1
  # Do community detection
  true_cluster <- c(rep(1, 29), rep(2, 21), rep(3, 30), rep(4, 20), 1, 3, 4)
  if (community_index == "AP"){
    res <- APC_fun(beta_mat)
  } else if(community_index == "LPA"){
    res <- LPA_fun(beta_mat, K = 5)
  } else if(community_index == "SCORE"){
    res <- SCORE_fun(beta_adj, K = 4)
  } else if (community_index == "SLIM"){
    res <- SLIM_fun(beta_mat_abs, K = 4)
  } else if (community_index == "M_SCORE"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SCORE", c_method = "kmeans")
  } else if (community_index == "M_SLIM"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SLIM", c_method = "kmeans")
  } else if (community_index == "M_SCORE_AP"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SCORE", c_method = "AP")
  } else if (community_index == "M_SCORE_LPA"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SCORE", c_method = "LPA")
  } else if (community_index == "M_SLIM_AP") {
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SLIM", c_method = "AP")
  } else if (community_index == "M_SLIM_LPA"){
    res <- MASE_fun(network_list, di = 4, K = 4, d_method = "SLIM", c_method = "LPA")
  }
  # return results
  return(fossil::rand.index(res, true_cluster))
}

## cor method
randindex_cor <- rep(NA, length(community_index))
for (i in 1:length(community_index)){
  randindex_cor[i] <- network_fun(dta_list, method_index = "cor", community_index = community_index[i])
}
saveRDS(randindex_cor, "out_cor_withnorm1000.rds")

cl <- makeCluster(20) # 4 workers
registerDoParallel(cl)

set.seed(0104)
out_res <- foreach(iter = 1:nrow(index)) %dorng% {
  try(network_fun(dta_list, method_index = "PCnet", npc_index = index[iter, 1], ncell_index = index[iter, 2], 
                  tensor_index = index[iter, 3], nq_index = index[iter, 4], community_index = index[iter, 5]))
}
stopCluster(cl)
saveRDS(out_res, "out_res_withnorm_1000.rds")
