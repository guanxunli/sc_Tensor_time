library(Matrix)
source("simulation_utility.R")

## load data
dta <- read.csv("SERGIO/diffOutput_S_1000.csv", header = FALSE)
rownames(dta) <- paste0("Gene", 1:nrow(dta))
colnames(dta) <- paste0("Cell", 1:ncol(dta))
nGenes <- nrow(dta)

dta_list <- list()
for (i in 1:10){
  dta_list[[i]] <- dta[, 1:(1000 * i)]
}

community_index <- c("AP", "SCORE", "LPA", "SLIM", "M_SCORE", "M_SLIM", "M_SCORE_AP", "M_SCORE_LPA", "M_SLIM_AP", "M_SLIM_LPA")

## define function
network_fun <- function(dta_list, community_index){
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
saveRDS(randindex_cor, "out_cor_withoutnorm.rds")

########################################## With new normalization method ##########################################
dta <- new_Normalization(dta)
dta_list <- list()
for (i in 1:10){
  dta_list[[i]] <- dta[, 1:(1000 * i)]
}

## cor method
randindex_cor <- rep(NA, length(community_index))
for (i in 1:length(community_index)){
  randindex_cor[i] <- network_fun(dta_list, method_index = "cor", community_index = community_index[i])
}
saveRDS(randindex_cor, "out_cor_withnorm.rds")

########################################################################################################
########################################## Check results ###############################################
########################################################################################################

library(ggplot2)
out_cor <- readRDS("results/correlation method/out_cor.rds")

# plot community detection method
png('results/community_method.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_cor) +
  geom_boxplot(aes(x = community_index, y = randindex, color = community_index)) +
  xlab("community_index") +
  ylab("rand index") +
  labs(title = "community method") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# plot whether normalization
png('results/normalization_method.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_cor) +
  geom_boxplot(aes(x = norm_index, y = randindex, color = norm_index)) +
  xlab("normalization or not") +
  ylab("rand index") +
  labs(title = "normalization method") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# plot how to generate dta
png('results/generate_method.png', width = 2500, height = 1500, res = 300)
ggplot(data = out_cor) +
  geom_boxplot(aes(x = dta, y = randindex, color = dta)) +
  xlab("generate method") +
  ylab("rand index") +
  labs(title = "generate method") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##################### Check the best method ################################
out_cor <- out_cor[order(out_cor$randindex, decreasing = TRUE), ]
head(out_cor)

library(ComplexHeatmap)
source("simulation_utility.R")

## load data
dta <- read.csv("SERGIO/diffOutput_S_1000.csv", header = FALSE)
rownames(dta) <- paste0("Gene", 1:nrow(dta))
colnames(dta) <- paste0("Cell", 1:ncol(dta))
nGenes <- nrow(dta)
dta <- new_Normalization(dta)
dta_list <- list()
for (i in 1:10){
  dta_list[[i]] <- dta[, (1000 * (i-1) + 1):(1000 * i)]
}
## make networks
network_list <- list()
for (i in 1:length(dta_list)){
  network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
  rownames(network_list[[i]]) <- NULL
  colnames(network_list[[i]]) <- NULL
}
Heatmap(network_list[[10]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

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
Heatmap(beta_mat, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

