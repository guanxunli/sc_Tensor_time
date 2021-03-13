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
res <- list()

#### No trajectary used
dta_list <- list()
time_vec <- 1:9
set.seed(1)
for (i in 1:9) {
  sample_index <- sample(seq_len(n_cell), 1000, replace = TRUE)
  dta_list[[i]] <- as.matrix(dta[, sample_index])
}

## make networks
res_notraj <- scTenifoldTime_beta(dta_list, time_vec, method = "cor", K = 5)
network_list <- list()
for (i in seq_len(length(dta_list))) {
  network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
  network_list[[i]] <- round(network_list[[i]], 2)
  print(paste0("Finish network", i))
}
set.seed(1)
network_tensor <- tensorDecomposition(network_list, K = 5, maxIter = 10000, maxError = 1e-5)
res$network_tensor_notraj <- network_tensor
res_regression <- my_regression(network_list = network_tensor, time_vec = time_vec)
beta_mat <- res_regression$beta_mat
rownames(beta_mat) <- rownames(dta)
set.seed(1)
E <- UMAP_order(dta = beta_mat)
res$E <- list()
res$E$no_traj <- E
print("Finish no trajectary part.")

#### Two basic method
## Do UMAP directly
set.seed(1)
E <- UMAP_order(dta = dta)
res$E$UMAP_dir <- E
## One network
dta_net <- round(cor(t(as.matrix(dta)), method = "spearman"), 2)
rownames(dta_net) <- rownames(dta)
set.seed(1)
E <- UMAP_order(dta = dta_net)
res$E$dta_net <- E
print("Finish two basic methods.")
saveRDS(res, "results_10000/res_cor1.rds")

#### Now consider our own method
dta_list <- list()
time_vec <- rep(NA, 9)
for (i in seq_len(8)) {
  dta_list[[i]] <- as.matrix(dta[, (500 * (i - 1) + 1):(500 * (i + 1))])
  time_vec[i] <- mean(dta_sudotime1[colnames(dta_list[[i]]), 3])
}
dta_list[[9]] <- as.matrix(dta[, 4001:n_cell])
time_vec[9] <- mean(dta_sudotime1[colnames(dta_list[[9]]), 3])
time_vec <- time_vec / max(time_vec)

## scTenifoldTime mehod
res_cor <- scTenifoldTime_beta(dta_list, time_vec, method = "cor", K = 10)
res$res_cor <- res_cor

## Without tensor decomposition
set.seed(1)
E <- UMAP_order(dta = res_cor$beta_mat)
res$E$beta_time <- E
beta_filter <- matrix(0, nrow = nrow(res_cor$beta_mat), ncol = nrow(res_cor$beta_mat))       
index <- which(res_cor$t_mat < 0.05)     
beta_filter[index] <- res_cor$beta_mat[index]  
rownames(beta_filter) <- colnames(beta_filter) <- rownames(dta_list[[1]])     
set.seed(1) 
E <- UMAP_order(dta = beta_filter)  
res$E$beta_time_filter <- E

# #### community detecyion
# res_APC <- APC_fun(beta_mat)
# res_SCORE <- SCORE_fun(beta_adj, K = length(unique(res_APC)))
# res_SLIM_AP <- MASE_fun(network_list, di = 50, K = 50, gamma = 0.25, m = 8)
# res$APC <- res_APC
# res$SCORE <- res_SCORE
# res$SLIM_AP <- res_SLIM_AP
# saveRDS(res, "results/res_cor1.rds")
# ## Transfer format
# gene_list <- rownames(dta)
# names(res_APC) <- gene_list
# names(res_SCORE) <- gene_list
# names(res_SLIM_AP) <- gene_list

# AP_list <- list()
# i <- 1
# tmp_community <- unique(res_APC)
# for (iter in 1:length(unique(res_APC))) {
#   index <- which(res_APC == tmp_community[iter])
#   if (length(index) > 5){
#     tmp <- names(res_APC)[index]
#     t_tmp <- as.numeric(t_mat[index, index])
#     if (min(t_tmp) < 0.05){
#       AP_list[[i]] <- tmp
#       i <- i + 1
#     }
#   }
# }

# SCORE_list <- list()
# i <- 1
# tmp_community <- unique(res_SCORE)
# for (iter in 1:length(unique(res_SCORE))) {
#   index <- which(res_SCORE == tmp_community[iter])
#   if (length(index) > 5){
#     tmp <- names(res_SCORE)[index]
#     t_tmp <- as.numeric(t_mat[index, index])
#     if (min(t_tmp) < 0.05){
#       SCORE_list[[i]] <- tmp
#       i <- i + 1
#     }
#   }
# }

# SLIM_AP_list <- list()
# i <- 1
# tmp_community <- unique(res_SLIM_AP)
# for (iter in 1:length(unique(res_SLIM_AP))) {
#   index <- which(res_SLIM_AP == tmp_community[iter])
#   if (length(index) > 5){
#     tmp <- names(res_SLIM_AP)[index]
#     t_tmp <- as.numeric(t_mat[index, index])
#     if (min(t_tmp) < 0.05){
#       SLIM_AP_list[[i]] <- tmp
#       i <- i + 1
#     }
#   }
# }

# res$APC_list <- AP_list
# res$SCORE_list <- SCORE_list
# res$SLIM_AP_list <- SLIM_AP_list
# saveRDS(res, "results/res_cor1.rds")

#### With tensor decomposition
set.seed(1)
E <- UMAP_order(dta = res_cor$beta_mat_tensor)
res$E$beta_time_tensor <- E
beta_filter <- matrix(0, nrow = nrow(res_cor$beta_mat_tensor), ncol = nrow(res_cor$beta_mat_tensor))       
index <- which(res_cor$t_mat_tensor < 0.05)     
beta_filter[index] <- res_cor$beta_mat_tensor[index]
rownames(beta_filter) <- colnames(beta_filter) <- rownames(dta_list[[1]])     
set.seed(1) 
E <- UMAP_order(dta = beta_filter)  
res$E$beta_time_tensor_filter <- E

## save results.
saveRDS(res, "results_7000/res_cor1.rds")
# saveRDS(res, "results_10000/res_cor1.rds")

# #### community detection
# res_APC <- APC_fun(beta_mat)
# res_SCORE <- SCORE_fun(beta_adj, K = length(unique(res_APC)))
# res_SLIM_AP <- MASE_fun(network_tensor, di = 50, K = 50, gamma = 0.25, m = 8)
# res$APC_tensor <- res_APC
# res$SCORE_tensor <- res_SCORE
# res$SLIM_AP_tensor <- res_SLIM_AP
# saveRDS(res, "results/res_cor1.rds")
# ## Transfer format
# gene_list <- rownames(dta)
# names(res_APC) <- gene_list
# names(res_SCORE) <- gene_list
# names(res_SLIM_AP) <- gene_list

# AP_list <- list()
# i <- 1
# tmp_community <- unique(res_APC)
# for (iter in 1:length(unique(res_APC))) {
#   index <- which(res_APC == tmp_community[iter])
#   if (length(index) > 5){
#     tmp <- names(res_APC)[index]
#     t_tmp <- as.numeric(t_mat[index, index])
#     if (min(t_tmp) < 0.05){
#       AP_list[[i]] <- tmp
#       i <- i + 1
#     }
#   }
# }

# SCORE_list <- list()
# i <- 1
# tmp_community <- unique(res_SCORE)
# for (iter in 1:length(unique(res_SCORE))) {
#   index <- which(res_SCORE == tmp_community[iter])
#   if (length(index) > 5){
#     tmp <- names(res_SCORE)[index]
#     t_tmp <- as.numeric(t_mat[index, index])
#     if (min(t_tmp) < 0.05){
#       SCORE_list[[i]] <- tmp
#       i <- i + 1
#     }
#   }
# }

# SLIM_AP_list <- list()
# i <- 1
# tmp_community <- unique(res_SLIM_AP)
# for (iter in 1:length(unique(res_SLIM_AP))) {
#   index <- which(res_SLIM_AP == tmp_community[iter])
#   if (length(index) > 5){
#     tmp <- names(res_SLIM_AP)[index]
#     t_tmp <- as.numeric(t_mat[index, index])
#     if (min(t_tmp) < 0.05){
#       SLIM_AP_list[[i]] <- tmp
#       i <- i + 1
#     }
#   }
# }

# res$APC_list_tensor <- AP_list
# res$SCORE_list_tensor <- SCORE_list
# res$SLIM_AP_list_tensor <- SLIM_AP_list
# saveRDS(res, "results/res_cor1.rds")


# ## check the final results.
# library(enrichR)
# library(igraph)

# score_function <- function(res_list, fdrThreshold = 0.05){
#   n <- length(res_list)
#   n_community <- 0
#   for (i in 1:length(res_list)){
#     gList <- res_list[[i]]
#     gList <- gList[!grepl('^mt-|^Rpl|^Rps',gList, ignore.case = TRUE)]
#     if (length(gList) == 0){
#       n_community <- n_community
#     } else{
#       E <- enrichr(gList, c('KEGG_2019_Mouse','Reactome_2016','WikiPathways_2019_Mouse'))
#       E <- do.call(rbind.data.frame, E)
#       E <- E[E$Adjusted.P.value < fdrThreshold,]
#       print(i)
#       if(isTRUE(nrow(E) > 0)){
#         n_community <- n_community  + 1
#       }
#     }
#   }
#   return(n_community / n)
# }

# res$APC_score <- score_function(res$APC_list)
# res$SCORE_score <- score_function(res$SCORE_list)
# res$SLIM_AP_socre <- score_function(res$SLIM_AP_list)
# saveRDS(res, "res_cor1.rds")
