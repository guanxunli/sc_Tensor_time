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
######### Community detection method ###################
########################################################

#### Affinity Propagation Clustering ####
APC_fun <- function(net){
  res_ap <- apcluster(s = negDistMat(), x = net)
  com_det_ap <- res_ap@clusters
  res <- rep(NA, nrow(net))
  for (i in 1:length(com_det_ap)){
    res[as.numeric(com_det_ap[[i]])] <- i
  }
  return(res)
}

#### SCORE method  ######################
# net is the network
# K is the number of communties
# Tn is the threshold
SCORE_fun <- function(net, K, Tn = NULL){
  # make net work symmetry
  net <- (net + t(net)) / 2
  net <- as(net, "sparseMatrix")
  # get the ratio matrix
  library(RSpectra)
  eig_net <- eigs_sym(net, k = K, which = "LM")$vectors
  eig_ratio <- matrix(eig_net[, 1], nrow = nrow(eig_net), ncol = K - 1)
  eig_ratio <- eig_net[, 2:K] / eig_ratio
  if (!is.null(Tn)){
    eig_ratio[which(eig_ratio > Tn)] <- Tn
    eig_ratio[which(eig_ratio < -Tn)] <- -Tn
  }
  # do k-means
  SCORE_kmeans <- kmeans(eig_ratio, centers = K, nstart = 10, iter.max = 50)
  return(SCORE_kmeans$cluster)
}

#### label propagation algorithm (LPA) ##
# net is the network
# K is the number of neighborhood
getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

LPA_fun <- function(net, K, itermax = 10000){
  # calculate KNN
  n <- nrow(net)
  net_knn <- get.knn(net, k = K)$nn.index
  res_label <- 1:n
  # label propagation
  for (iter in 1:itermax) {
    res_label_new <- rep(NA, n)
    for (i in 1:n){
      res_label_new[i] <- as.numeric(getmode(res_label[c(net_knn[i, ], i)]))
    }
    if (sum(abs(res_label_new - res_label)) == 0) break
    
    res_label <- res_label_new
  }
  res <- list()
  res$iter <- iter
  res$community <- res_label
  return(res)
}

#### non-negative matrix factorlization (NMF) ##
# net is the network
# K is the number of communities
NMF_fun <- function(net, K){
  net <- abs(net)
  res_NMF <- nmf(x = net, rank = K)@fit@W
  res_NMF <- apply(res_NMF, 1, which.max)
  return(res_NMF)
}

#### Symmetrized Laplacian Inverse Matrix (SLIM)
# net is the network
# K is the number of communities
# m is the number of terms of summation
# gamma is decrease rate
SLIM_fun <- function(net, gamma = 0.25, m = 8, K = 3){
  # Calculate the inverse Laplacian matrix
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
  
  # Spectral decomposition
  M_eig <- RSpectra::eigs_sym(M_hat, k = K, which = "LM")
  X_hat <- M_eig$vectors
  
  # do k-means
  SLIM_kmeans <- kmeans(X_hat, centers = K, nstart = 10, iter.max = 50)
  return(SLIM_kmeans$cluster)
}

#### Multiple adjacency spectral embedding (MASE)
MASE_fun <- function(net_list, di = 5, K = di){
  # Do adjacency spectral embedding
  m <- length(net_list)
  U <- matrix(NA, nrow = nrow(net_list[[1]]), ncol = m * di)
  for (i in 1:m){
    A <- net_list[[i]]
    A_eig <- RSpectra::eigs_sym(A, k = di, which = "LM")
    U[ ,(1 + di * (i - 1)):(di * i)] <- A_eig$vectors
  }
  U_svd <- RSpectra::svds(U, k = K)
  V <- U_svd$u
  
  # Generate R
  Rhat_list <- list()
  for (i in 1:m){
    Rhat_list[[i]] <- crossprod(V, net_list[[i]] %*% V)
  }
  
  # return results
  res <- list()
  res$V <- V
  res$R <- Rhat_list
  return(res)
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