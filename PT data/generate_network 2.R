library(Matrix)
source("utility.R")

## load data
dta <- readRDS("data/dta.rds")
dta_sudotime <- read.csv("data/monocle_rank.csv")
dta_sudotime <- dta_sudotime[order(dta_sudotime$CCAT_score), ]
dta <- dta[, dta_sudotime$Barcode]

## previous method
dta_net_old <- scTeni_makeNet(dta, norm_method = "old", nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                          nc_q = 0.05, tensor_dc = "FALSE", td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5)
saveRDS(dta_net_old, "dta_net_old.rds")
dta_net_old_td <- scTeni_makeNet(dta, norm_method = "old", nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                              nc_q = 0.05, tensor_dc = "TRUE", td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5)
saveRDS(dta_net_old_td, "dta_net_old_td.rds")
dta_net_new <- scTeni_makeNet(dta, norm_method = "new", nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                              nc_q = 0.05, tensor_dc = "FALSE", td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5)
saveRDS(dta_net_new, "dta_net_new.rds")
dta_net_new_td <- scTeni_makeNet(dta, norm_method = "new", nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                              nc_q = 0.05, tensor_dc = "TRUE", td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5)
saveRDS(dta_net_new_td, "dta_net_new_td.rds")

## 8 data set
dta_list <- list()
for (i in 1:8){
  if (i == 8){
    dta_list[[i]] <- dta[, (543 * (i - 1) + 1) : ncol(dta)]
  } else{
    dta_list[[i]] <- dta[, (543 * (i - 1) + 1) : (543 * i)]
  }
}

## generate old networks
Network_list_old <- list()
for (i in 1:8){
  Network_list_old[[i]] <- scTeni_makeNet(scRNAsq_list[[i]], norm_method = "old", tensor_dc = "FALSE")
}
nGenes <- nrow(Network_list_old[[1]])
for(i in 1:8){
  png(paste0('oMatrix_old',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_old[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_old)
saveRDS(reg_res, "reg_res_old.rds")

Network_list_old_td <- list()
for (i in 1:8){
  Network_list_old_td[[i]] <- scTeni_makeNet(scRNAsq_list[[i]],norm_method = "old", tensor_dc = "TRUE", td_K = 20)
}

for(i in 1:8){
  png(paste0('tMatrix_old',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_old_td[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_old_td)
saveRDS(reg_res, "reg_res_old_td.rds")

## generate new networks
Network_list_new <- list()
for (i in 1:8){
  Network_list_new[[i]] <- scTeni_makeNet(scRNAsq_list[[i]], norm_method = "new", tensor_dc = "FALSE")
}

for(i in 1:8){
  png(paste0('oMatrix_new',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_new[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_new)
saveRDS(reg_res, "reg_res_new.rds")

Network_list_new_td <- list()
for (i in 1:8){
  Network_list_new_td[[i]] <- scTeni_makeNet(scRNAsq_list[[i]],norm_method = "new", tensor_dc = "TRUE",td_K = 20)
}

for(i in 1:8){
  png(paste0('tMatrix_new',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_new_td[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_new_td)
saveRDS(reg_res, "reg_res_new_td.rds")
