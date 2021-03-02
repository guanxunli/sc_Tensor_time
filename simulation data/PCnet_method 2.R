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

########################################################################################################
##################### generate networks ################################################################
########################################################################################################

npc_index <- c(3, 5, 9, 13)
ncell_index <- c(500, 800, 1000, 1200, 1500)
nq_index <- c(0.05, 0.25, 0.5)
tensor_index <- c(0, 3, 20)
index <- expand.grid(npc_index = npc_index, ncell_index = ncell_index,
                     tensor_index = tensor_index, nq_index = nq_index)

## define function
network_fun <- function(dta_list, npc_index, ncell_index, tensor_index, nq_index){
  ## make networks
  network_list <- list()
  for (i in 1:length(dta_list)){
    network_list[[i]] <- scTeni_makeNet(dta_list[[i]], norm_method = "old", nc_nNet = 10, nc_nCells = ncell_index, nc_nComp = npc_index, 
                                        nc_q = nq_index, td_K = tensor_index)
  }
  return(network_list)
}

## Do simulation
library(foreach)
library(doParallel)
library(doRNG)

cl <- makeCluster(20) # 4 workers
registerDoParallel(cl)

set.seed(0104)
out_res <- foreach(iter = 1:nrow(index)) %dorng% {
  network_fun(dta_list, npc_index = index[iter, 1], ncell_index = index[iter, 2], 
              tensor_index = index[iter, 3], nq_index = index[iter, 4])
}
stopCluster(cl)
saveRDS(out_res, "out_res_withoutnorm.rds")

#### With new normalization method 
dta <- new_Normalization(dta)
dta_list <- list()
for (i in 1:10){
  dta_list[[i]] <- dta[, (1000 * (i - 1) + 1):(1000 * i)]
}

## define function
network_fun <- function(dta_list, npc_index, ncell_index, tensor_index, nq_index){
  ## make networks
  network_list <- list()
  for (i in 1:length(dta_list)){
    network_list[[i]] <- scTeni_makeNet(dta_list[[i]], norm_method = "new", nc_nNet = 10, nc_nCells = ncell_index, nc_nComp = npc_index, 
                                        nc_q = nq_index, td_K = tensor_index)
  }
  return(network_list)
}

cl <- makeCluster(20) # 4 workers
registerDoParallel(cl)

set.seed(0104)
out_res <- foreach(iter = 1:nrow(index)) %dorng% {
  network_fun(dta_list, npc_index = index[iter, 1], ncell_index = index[iter, 2], 
              tensor_index = index[iter, 3], nq_index = index[iter, 4])
}
stopCluster(cl)
saveRDS(out_res, "out_res_withnorm.rds")

########################################################################################################
################### community detection ################################################################
########################################################################################################
library(Matrix)
source("simulation_utility.R")

## community detection
community_detect <- function(net_list){
  # Getting beta mat
  n_net <- length(net_list)
  X_mat <- cbind(1, 1:n_net)
  tmp <- as.numeric(net_list[[1]])
  Y <- tmp
  for (i in 2:10){
    tmp <- as.numeric(net_list[[i]])
    Y <- rbind(Y, tmp)
  }
  ## Do regress
  beta_coef <- solve(crossprod(X_mat), crossprod(X_mat, Y))
  beta_mat <- matrix(beta_coef[2, ], nrow = nrow(net_list[[1]]))
  beta_mat_abs <- abs(beta_mat)
  beta_adj <- round(beta_mat_abs, 2)
  beta_adj[beta_adj > 0] <- 1
  # Do community detection
  true_cluster <- c(rep(1, 29), rep(2, 21), rep(3, 30), rep(4, 20), 1, 3, 4)
  rand_index <- list()
  res <- APC_fun(beta_mat)
  rand_index$AP <- fossil::rand.index(res, true_cluster)
  res <- LPA_fun(beta_mat, K = 5)
  rand_index$LPA <- fossil::rand.index(res, true_cluster)
  res <- SCORE_fun(beta_adj, K = 4)
  rand_index$SCORE <- fossil::rand.index(res, true_cluster)
  res <- SLIM_fun(beta_mat_abs, K = 4)
  rand_index$SLIM <- fossil::rand.index(res, true_cluster)
  res <- MASE_fun(net_list, di = 5, K = 4)
  rand_index$SCORE_kmeans <- fossil::rand.index(res$SCORE_kmeans, true_cluster) 
  rand_index$SCORE_APC <- fossil::rand.index(res$SCORE_APC, true_cluster) 
  rand_index$SCORE_LPA <- fossil::rand.index(res$SCORE_LPA, true_cluster) 
  rand_index$SLIM_kmeans <- fossil::rand.index(res$SLIM_kmeans, true_cluster) 
  rand_index$SLIM_APC <- fossil::rand.index(res$SLIM_APC, true_cluster) 
  rand_index$SLIM_LPA <- fossil::rand.index(res$SLIM_LPA, true_cluster) 
  # return
  return(rand_index)
}

library(foreach)
library(doParallel)
library(doRNG)

## with normalization
net_res <- readRDS("out_res_withnorm_1000.rds")

cl <- makeCluster(20) # 4 workers
registerDoParallel(cl)

set.seed(0104)
out_res <- foreach(iter = 1:length(net_res)) %dorng% {
  community_detect(net_res[[iter]])
}
stopCluster(cl)
saveRDS(out_res, "out_res_withnorm_1000.rds")

net_res <- readRDS("out_res_withoutnorm_1000.rds")
## without normalization
cl <- makeCluster(20) # 4 workers
registerDoParallel(cl)

set.seed(0104)
out_res <- foreach(iter = 1:length(net_res)) %dorng% {
  community_detect(net_res[[iter]])
}
stopCluster(cl)
saveRDS(out_res, "out_res_withoutnorm_1000.rds")

########################################################################################################
########################################## Check results ###############################################
########################################################################################################
# out_res <- readRDS("results/PCnet method/out_res_withnorm_1000.rds")
# npc_index <- as.factor(c(3, 5, 9, 13))
# ncell_index <- as.factor(c(500, 700, 800, 1000))
# nq_index <- as.factor(c(0.05, 0.25, 0.5))
# tensor_index <- as.factor(c(0, 3, 20))
# community_method <- names(out_res[[1]])
# index <- expand.grid(community_method = community_method, npc_index = npc_index, ncell_index = ncell_index,
#                      tensor_index = tensor_index, nq_index = nq_index)
# res_df_withnorm_1000 <- as.data.frame(index)
# res_df_withnorm_1000$rand_index <- as.numeric(unlist(out_res))
# 
# out_res <- readRDS("results/PCnet method/out_res_withnorm_2000.rds")
# ncell_index <- as.factor(c(500, 800, 1000, 1200, 1500))
# index <- expand.grid(community_method = community_method, npc_index = npc_index, ncell_index = ncell_index,
#                      tensor_index = tensor_index, nq_index = nq_index)
# res_df_withnorm_2000 <- as.data.frame(index)
# res_df_withnorm_2000$rand_index <- as.numeric(unlist(out_res))
# 
# out_res <- readRDS("results/PCnet method/out_res_withnorm.rds")
# index <- expand.grid(community_method = community_method, npc_index = npc_index, ncell_index = ncell_index,
#                      tensor_index = tensor_index, nq_index = nq_index)
# res_df_withnorm <- as.data.frame(index)
# res_df_withnorm$rand_index <- as.numeric(unlist(out_res))
# 
# out_res <- readRDS("results/PCnet method/out_res_withoutnorm_1000.rds")
# ncell_index <- as.factor(c(500, 700, 800, 1000))
# index <- expand.grid(community_method = community_method, npc_index = npc_index, ncell_index = ncell_index,
#                      tensor_index = tensor_index, nq_index = nq_index)
# res_df_withoutnorm_1000 <- as.data.frame(index)
# res_df_withoutnorm_1000$rand_index <- as.numeric(unlist(out_res))
# 
# out_res <- readRDS("results/PCnet method/out_res_withoutnorm_2000.rds")
# ncell_index <- as.factor(c(500, 800, 1000, 1200, 1500))
# index <- expand.grid(community_method = community_method, npc_index = npc_index, ncell_index = ncell_index,
#                      tensor_index = tensor_index, nq_index = nq_index)
# res_df_withoutnorm_2000 <- as.data.frame(index)
# res_df_withoutnorm_2000$rand_index <- as.numeric(unlist(out_res))
# 
# out_res <- readRDS("results/PCnet method/out_res_withoutnorm.rds")
# index <- expand.grid(community_method = community_method, npc_index = npc_index, ncell_index = ncell_index,
#                      tensor_index = tensor_index, nq_index = nq_index)
# res_df_withoutnorm <- as.data.frame(index)
# res_df_withoutnorm$rand_index <- as.numeric(unlist(out_res))
# 
# res_df <- rbind(res_df_withnorm_1000, res_df_withnorm_2000, res_df_withnorm,
#                 res_df_withoutnorm_1000, res_df_withoutnorm_2000, res_df_withoutnorm)
# res_df$norm_index <- rep(c("TRUE", "FALSE"), each = 5040)
# res_df$dta <- rep(as.factor(c(rep(1000, 1440), rep(2000, 1800), rep(0, 1800))), 2)

# saveRDS(res_df, "results/PCnet method/out_PCnet.rds")
out_PCnet <- readRDS("results/PCnet method/out_PCnet.rds")
library(ggplot2)

png('results/PCnet method/community_method.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_PCnet) +
  geom_boxplot(aes(x = community_method, y = rand_index, color = community_method)) +
  xlab("community_index") +
  ylab("rand index") +
  labs(title = "community method") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('results/PCnet method/number_PC.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_PCnet) +
  geom_boxplot(aes(x = npc_index, y = rand_index, color = npc_index)) +
  xlab("npc_index") +
  ylab("rand index") +
  labs(title = "number of PCs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('results/PCnet method/number_cells.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_PCnet) +
  geom_boxplot(aes(x = ncell_index, y = rand_index, color = ncell_index)) +
  xlab("ncell_index") +
  ylab("rand index") +
  labs(title = "number of cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('results/PCnet method/tensor_index.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_PCnet) +
  geom_boxplot(aes(x = tensor_index, y = rand_index, color = tensor_index)) +
  xlab("tensor_index") +
  ylab("rand index") +
  labs(title = "tensor influence") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('results/PCnet method/quantile_index.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_PCnet) +
  geom_boxplot(aes(x = nq_index, y = rand_index, color = nq_index)) +
  xlab("nq_index") +
  ylab("rand index") +
  labs(title = "quantile influence") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('results/PCnet method/normalization_method.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_PCnet) +
  geom_boxplot(aes(x = norm_index, y = rand_index, color = norm_index)) +
  xlab("norm_index") +
  ylab("rand index") +
  labs(title = "normalization influence") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png('results/PCnet method/generate_method.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_PCnet) +
  geom_boxplot(aes(x = dta, y = rand_index, color = dta)) +
  xlab("group method") +
  ylab("rand index") +
  labs(title = "group method") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##################### Check the best method ################################
out_PCnet <- out_PCnet[order(out_PCnet$rand_index, decreasing = TRUE), ]
head(out_PCnet, 20)

#### check the best regression method
library(Matrix)
library(ComplexHeatmap)
source("simulation_utility.R")
dta <- read.csv("SERGIO/diffOutput_S_1000.csv", header = FALSE)
rownames(dta) <- paste0("Gene", 1:nrow(dta))
colnames(dta) <- paste0("Cell", 1:ncol(dta))
nGenes <- nrow(dta)

dta <- new_Normalization(dta)
network_list <- list()
for (i in 1:10){
  dta_tmp <- dta[, (1000 * (i - 1) + 1):(1000 * i)]
  network_list[[i]] <- scTeni_makeNet(dta_tmp, norm_method = "new", nc_nNet = 10, nc_nCells = 1000, nc_nComp = 5,
                                 nc_q = 0.5, td_K = 3)
  rownames(network_list[[i]]) <- NULL
  colnames(network_list[[i]]) <- NULL
}
Heatmap(network_list[[10]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
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
Heatmap(beta_mat, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

#### Check the best SLIM method
network_list <- list()
for (i in 1:10){
  dta_tmp <- dta[, 1:(1000 * i)]
  network_list[[i]] <- scTeni_makeNet(dta_tmp, norm_method = "new", nc_nNet = 10, nc_nCells = 1200, nc_nComp = 5,
                                      nc_q = 0.5, td_K = 3)
  rownames(network_list[[i]]) <- NULL
  colnames(network_list[[i]]) <- NULL
}
Heatmap(network_list[[10]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

## Do SLIM transformation
network_SLIM <- lapply(network_list,  function(net, gamma = 0.25, m = 8, K = 3){
  # Calculate the inverse Laplacian matrix
  net <- abs(net)
  D <- rowSums(net)
  net[which(D == 0), ] <- 1/ncol(net)
  D <- rowSums(net)
  DinvA <- net / D
  W_hat <- 0
  alpha <- exp(-gamma)
  tmp_c <- alpha
  tmp_m <- DinvA
  for (i in 1:m){
    W_hat <- W_hat + tmp_c * tmp_m
    tmp_c <- tmp_c * tmp_c
    tmp_m <- tmp_m %*% tmp_m
  }
  
  # Make it symmetry
  M_hat <- (W_hat + t(W_hat)) / 2
  diag(M_hat) <- 0
  return(M_hat)
})
Heatmap(network_SLIM[[8]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

