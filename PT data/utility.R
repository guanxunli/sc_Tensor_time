###########################################################
######## Regression function ##############################
###########################################################

## define regression function
my_regression <- function(X, Y){
  n <- dim(X)[2]
  df <- dim(X)[1] - n
  # Do regress
  beta_coef <- solve(crossprod(X), crossprod(X, Y))
  
  # t-test
  lm_t_p <- matrix(NA, nrow = n - 1, ncol = dim(Y)[2])
  index <- which(beta_coef[1, ] != 0)
  Y_fit <- X %*% beta_coef[, index]
  sigma_fit <- colSums((Y[, index] - Y_fit)^2) / df
  
  
  for (i in 2:n){
    sd_t <- sqrt(solve(crossprod(X))[i, i] * sigma_fit)
    t_test <- beta_coef[i, index] / sd_t
    t_test_p <- 2 * (1 - pt(abs(t_test), df = df))
    lm_t_p[i - 1, index] <- t_test_p
  }
  
  # # F-test
  # X_F <- X[, -ncol(X), drop=FALSE]
  # beta_coef_F <- solve(crossprod(X_F), crossprod(X_F, Y))
  # Y_fit_F <- X_F %*% beta_coef_F[, index]
  # sigma_fit_F <- colSums((Y[, index] - Y_fit_F)^2) / df
  # F_test <- (sigma_fit_F - sigma_fit) * df / sigma_fit 
  # F_test_p <- 1 - pf(F_test, df1 = 1, df2 = df)
  # lm_F_p <- rep(NA, ncol = dim(Y)[2])
  # lm_F_p[index] <- F_test_p
  
  res <- list()
  res$coef <- beta_coef
  res$lm_t_p <- lm_t_p
  # res$lm_F_p <- lm_F_p
  return(res)
}

## Input network list
regression_res <- function(network_list, method = "linear"){
  n_list <- length(network_list)
  mat_Y <- as.numeric(network_list[[1]])
  for (i in 2:n_list){
    mat_Y <- rbind(mat_Y, as.numeric(network_list[[i]]))
  }
  mat_Y <- as.matrix(mat_Y)
  
  if (method == "linear"){
    X <- matrix(c(rep(1, 10), (1:10)), nrow = 10)
  } else{
    X <- matrix(c(rep(1, 10), (1:10), (1:10)^2), nrow = 10)
  }
  lm_res <- my_regression(X, mat_Y)
  return(lm_res)
}

########################################################
#### scTenifoldNet function modification. ##############
########################################################

new_Normalization <- function(X){
  X_sum <- sum(X)
  X_rowsum <- rowSums(X)
  X_colSum <- colSums(X)
  u_hat <- X_rowsum %*% t(X_colSum) / X_sum
  Z <- (X - u_hat) / sqrt(u_hat + u_hat^2 / 100)
  return(Z)
}

list2network <- function(x) {
  n_list <- length(x)
  x_net <- 0
  for (i in n_list){
    x_net <- x_net + x[[i]]
  }
  x_net <- as.matrix(x_net / n_list)
  rownames(x_net) <- NULL
  colnames(x_net) <- NULL
  return(x_net)
}

## make networks
scTeni_makeNet <- function(X, norm_method = "new", nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                          nc_q = 0.05, tensor_dc = "TRUE", td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5){
  
  X <- scTenifoldNet::scQC(X, minLibSize = 1000, removeOutlierCells = TRUE, minPCT = 0.05, maxMTratio = 0.1)
  
  # Counts per million (CPM) normalization
  if (norm_method == "new"){
    X <- new_Normalization(X)
  } else{
    X <- scTenifoldNet::cpmNormalization(X)
  }

  # Comparing gene ids.
  xNames <- rownames(X)
  nGenes <- length(xNames)
  
  # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
  set.seed(1)
  xList <- scTenifoldNet::makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
  
  # CANDECOMP/PARAFRAC Tensor Decomposition
  if (tensor_dc == TRUE){
    set.seed(1)
    tensorOut <- scTenifoldNet::tensorDecomposition(xList, K = td_K, maxIter = td_maxIter, maxError = td_maxError)
    
    # Split of tensor output
    tX <- as.matrix(tensorOut$X)
    rownames(tX) <- NULL
    colnames(tX) <- NULL
  } else{
    tX <- list2network(xList)
  }

  # Return
  return(tX)
}