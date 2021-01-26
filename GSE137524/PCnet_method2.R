library(Matrix)
source("utility.R")

#### Sudotime 2
## load data
dta <- readRDS("data/dta_raw.rds")
dta_sudotime2 <- read.csv("data/Ordered_SCC6_repeat2.csv")
dta_sudotime2 <- dta_sudotime2[order(dta_sudotime2$O2), ]
dta <- dta[, dta_sudotime2$X]
dta <- scTenifoldNet::scQC(as.matrix(dta))
dta <- new_Normalization(dta)

dta_list <- list()
for (i in 1:6){
  dta_list[[i]] <- as.matrix(dta[, (400 * (i - 1) + 1):(400 * (i + 1))])
}
dta_list[[7]] <- as.matrix(dta[, 2401:3241])

## make networks
network_list <- list()
for (i in 1:length(dta_list)){
  network_list[[i]] <- scTeni_makeNet(dta_list[[i]], norm_method = "new", nc_nNet = 10, nc_nCells = 500, nc_nComp = 5, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                                                  nc_q = 1, td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5)
}

res <- list()
res$network_list <- network_list
saveRDS(res, "res_PCnet2.rds")

## Getting beta mat
n_net <- length(network_list)
X_mat <- cbind(1, 1:n_net)
tmp <- as.numeric(network_list[[1]])
nGenes <- nrow(network_list[[1]])
Y <- tmp
for (i in 2:n_net){
  tmp <- as.numeric(network_list[[i]])
  Y <- rbind(Y, tmp)
}
## Do regress
beta_coef <- solve(crossprod(X_mat), crossprod(X_mat, Y))
beta_mat <- matrix(beta_coef[2, ], nrow = nGenes)
beta_adj <- round(abs(beta_mat), 2)
beta_adj[which(beta_adj > 0)] <- 1

## t-test
lm_t_p <- matrix(NA, nrow = ncol(X_mat) - 1, ncol = dim(Y)[2])
df <- nrow(X_mat) - ncol(X_mat)
index <- which(beta_coef[1, ] != 0)
Y_fit <- X_mat %*% beta_coef[, index]
sigma_fit <- colSums((Y[, index] - Y_fit)^2) / df

for (i in 2:ncol(X_mat)){
  sd_t <- sqrt(solve(crossprod(X_mat))[i, i] * sigma_fit)
  t_test <- beta_coef[i, index] / sd_t
  t_test_p <- 2 * (1 - pt(abs(t_test), df = df))
  lm_t_p[i - 1, index] <- t_test_p
}
t_mat <- matrix(lm_t_p, nrow = nGenes)

res <- list()
res$beta_mat <- beta_mat
res$t_mat <- t_mat
saveRDS(res, "res_PCnet2.rds")

## community detecyion
res_APC <- APC_fun(beta_mat)
res_SLIM_AP <- MASE_fun(network_list, di = 50, K = 50, gamma = 0.25, m = 8)

res$APC <- res_APC
res$SLIM_AP <- res_SLIM_AP
saveRDS(res, "res_cor1.rds")

gene_list <- rownames(dta)
names(res_APC) <- gene_list
names(res_SLIM_AP) <- gene_list

AP_list <- list()
i <- 1
tmp_community <- unique(res_APC)
for (iter in 1:length(unique(res_APC))) {
  index <- which(res_APC == tmp_community[iter])
  if (length(index) > 5){
    tmp <- names(res_APC)[index]
    t_tmp <- as.numeric(t_mat[index, index])
    t_tmp <- na.omit(t_tmp)
    if (-is.null(t_tmp)){
      if (min(t_tmp) < 0.05){
        AP_list[[i]] <- tmp
        i <- i + 1
      }
    }
  }
}

SLIM_AP_list <- list()
i <- 1
tmp_community <- unique(res_SLIM_AP)
for (iter in 1:length(unique(res_SLIM_AP))) {
  index <- which(res_SLIM_AP == tmp_community[iter])
  if (length(index) > 5){
    tmp <- names(res_SLIM_AP)[index]
    t_tmp <- as.numeric(t_mat[index, index])
    t_tmp <- na.omit(t_tmp)
    if (-is.null(t_tmp)){
      if (min(t_tmp) < 0.05){
        SLIM_AP_list[[i]] <- tmp
        i <- i + 1
      }
    }
  }
}

res$APC <- AP_list
res$SLIM_AP <- SLIM_AP_list
saveRDS(res, "res_PCnet2.rds")

## check the final results.
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
      print(i)
      if(isTRUE(nrow(E) > 0)){
        n_community <- n_community  + 1
      }
    }
  }
  return(n_community / n)
}

res_test <- list()
res_test[[1]] <- res$APC
res_test[[2]] <- res$SLIM_AP
res$score <- lapply(res_test, score_function)
saveRDS(res, "res_PCnet2.rds")
