library(Matrix)
library(ComplexHeatmap)
library(ggplot2)

dta <- read.csv("SERGIO/diffOutput_S_1000.csv", header = FALSE)
range(dta)
dim(dta)
rownames(dta) <- paste0("Gene", 1:nrow(dta))
colnames(dta) <- paste0("Cell", 1:ncol(dta))
source("utility.R")
nGenes <- nrow(dta)
#############################################################################
############## Using X_{t} to predict X_{t} #################################
#############################################################################

# png(paste0('network_',nComp,'q_', q, 'nCells_', nCells, '.png'), width = 1000, height = 1000, res = 300)
# print(Heatmap(dta_net, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
# dev.off()

#### All data sets
## old method
res <- scTeni_makeNet(dta, norm_method = "old", nc_nCells = 500, nc_nComp = 3, nc_q = 0.3, td_K = 3)
Heatmap(res$oX, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
Heatmap(res$tX, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

## new method
res <- scTeni_makeNet(dta, norm_method = "new", nc_nCells = 500, nc_nComp = 3, nc_q = 0.3, td_K = 3)
Heatmap(res$oX, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
Heatmap(res$tX, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

## quick test
res <- scTeni_makeNet(dta[,1:1000], norm_method = "old", nc_nCells = 500, nc_nComp = 3, nc_q = 0.3, td_K = 3)
Heatmap(res$oX %*% res$oX, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
Heatmap(res$tX, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

## correlation
res_cor <- cor(t(dta[, 5000:10000]), method = "spearman")
rownames(res_cor) <- NULL
colnames(res_cor) <- NULL
diag(res_cor) <- 0
# abs_res <- abs(res_cor)
# res_cor[abs_res < quantile(abs_res, 0.75)] <- 0
Heatmap(res_cor, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

#### 10 data set
## divide into different data sets
scRNAsq_list <- list()
for (i in 1:10){
  scRNAsq_list[[i]] <- dta[, (1 + 1000 * (i - 1)):(1000 * i)]
}

## generate old networks
Network_list_old <- list()
for (i in 1:10){
  Network_list_old[[i]] <- scTeni_makeNet(scRNAsq_list[[i]], norm_method = "old", nc_q = 0.45, tensor_dc = "FALSE", nc_nCells = 800, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('oMatrix_old',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_old[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_old)
saveRDS(reg_res, "reg_res_old.rds")

Network_list_old_td <- list()
for (i in 1:10){
  Network_list_old_td[[i]] <- scTeni_makeNet(scRNAsq_list[[i]],norm_method = "old", nc_q = 0.45, tensor_dc = "TRUE", nc_nCells = 800, td_K = 20, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('tMatrix_old',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_old_td[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_old_td)
saveRDS(reg_res, "reg_res_old_td.rds")

## generate new networks
Network_list_new <- list()
for (i in 1:10){
  Network_list_new[[i]] <- scTeni_makeNet(scRNAsq_list[[i]], norm_method = "new", nc_q = 0.45, tensor_dc = "FALSE", nc_nCells = 800, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('oMatrix_new',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_new[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_new)
saveRDS(reg_res, "reg_res_new.rds")

Network_list_new_td <- list()
for (i in 1:10){
  Network_list_new_td[[i]] <- scTeni_makeNet(scRNAsq_list[[i]],norm_method = "new", nc_q = 0.45, tensor_dc = "TRUE", nc_nCells = 800, td_K = 20, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('tMatrix_new',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_new_td[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_new_td)
saveRDS(reg_res, "reg_res_new_td.rds")

#### Community detection
plot_coef <- function(x_rds){
  reg_beta <- x_rds$coef
  reg_beta <- matrix(reg_beta[2, ], nrow = nGenes, ncol = nGenes)
  print(Heatmap(reg_beta, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  reg_beta_pca <- prcomp(reg_beta, scale. = FALSE)
  print(ggplot(data = NULL, aes(x = reg_beta_pca$x[, 1], y = reg_beta_pca$x[, 2])) +
    geom_point(color = as.factor(c(rep(1, 28), rep(2, 21), rep(3, 30), rep(4, 20), rep(5, 3)))) +
    labs(x = "PC1", y = "PC2", title = "PCA results"))
  reg_beta_tsne <- Rtsne::Rtsne(reg_beta)
  reg_beta_tsne <- reg_beta_tsne$Y
  print(ggplot(data = NULL, aes(x = reg_beta_tsne[, 1], y = reg_beta_tsne[, 2])) +
          geom_point(color = as.factor(c(rep(1, 28), rep(2, 21), rep(3, 30), rep(4, 20), rep(5, 3)))) +
          labs(x = "tsne1", y = "tsne2", title = "tsne results"))
}
reg_res <- readRDS("results/simulation/reg_res_old.rds")
plot_coef(reg_res)
reg_res <- readRDS("results/simulation/reg_res_old_td.rds")
plot_coef(reg_res)
reg_res <- readRDS("results/simulation/reg_res_new.rds")
plot_coef(reg_res)
reg_res <- readRDS("results/simulation/reg_res_new_td.rds")
plot_coef(reg_res)

