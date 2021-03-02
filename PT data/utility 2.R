########################################################
######### Community detection method ###################
########################################################

#### Affinity Propagation Clustering ####
APC_fun <- function(net){
  res_ap <- apcluster::apcluster(s = apcluster::negDistMat(), x = net)
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
  eig_net <- RSpectra::eigs_sym(net, k = K, which = "LM")$vectors
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


#### Multiple adjacency spectral embedding (MASE)
MASE_fun <- function(net_list, di = 10, K = 10, gamma = 0.25, m = 8){
  # Do dimensional reduction SLIM
  num_net <- length(net_list)
  U <- matrix(NA, nrow = nrow(net_list[[1]]), ncol = num_net * di)
  for (n_iter in 1:num_net){
    net <- net_list[[n_iter]]
    net <- abs(net)
    # Calculate the inverse Laplacian matrix
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
    # Spectral decomposition
    M_eig <- RSpectra::eigs_sym(M_hat, k = di, which = "LM")
    X_hat <- M_eig$vectors
    U[, (di * (n_iter - 1) + 1):(di * n_iter)] <- X_hat
  }
  U_svd <- RSpectra::svds(U, k = K)
  V <- U_svd$u
  res_SLIM_APC <- APC_fun(V)
  # return
  return(res_SLIM_APC)
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

pcNet <- function(X,
                  nComp = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0, verbose = TRUE) {
  xClass <- class(X)[[1]]
  validClass <- xClass %in% c('matrix', 'dgCMatrix')
  if (!validClass) {
    stop('Input should be a matrix with cells as columns and genes as rows')
  }
  if (nComp < 2 | nComp >= nrow(X)) {
    stop('nCom should be greater or equal than 2 and lower than the total number of genes')
  }
  gNames <- rownames(X)
  pcCoefficients <- function(K) {
    # Taking out the gene to be regressed out
    y <- X[, K]
    Xi <- X
    Xi <- Xi[, -K]
    # Step 1: Perform PCA on the observed covariates data matrix to obtain $n$ number of the principal components.
    coeff <- RSpectra::svds(Xi, nComp)$v
    score <- Xi %*% coeff
    # Step 2: Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients.
    score <-
      Matrix::t(Matrix::t(score) / (apply(score, 2, function(X) {
        sqrt(sum(X ^ 2))
      }) ^ 2))
    # Step 3: Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)
    
    return(Beta)
  }
  
  # # Standardizing the data
  # X <- (scale(Matrix::t(X)))
  X <- t(X)
  
  # Identify the number of rows in the input matrix
  n <- ncol(X)
  
  # Generate the output matrix
  A <- 1 - diag(n)
  
  # Apply the principal component regression for each gene
  if(verbose){
    B <- pbapply::pbsapply(seq_len(n), pcCoefficients)  
  } else {
    B <- sapply(seq_len(n), pcCoefficients)  
  }
  
  # Transposition of the Beta coefficient matrix
  B <- t(B)
  
  # Replacing the values in the output matrix
  for (K in seq_len(n)) {
    A[K, A[K, ] == 1] = B[K, ]
  }
  
  # Making the output matrix symmetric
  if (isTRUE(symmetric)) {
    A <- (A + t(A)) / 2
  }
  
  # Absolute values for scaling and filtering
  absA <- abs(A)
  
  # Scaling the output matrix
  if (isTRUE(scaleScores)) {
    A <- (A / max(absA))
  }
  
  # Filtering the output matrix
  A[absA < quantile(absA, q)] <- 0
  
  # Setting the diagonal to be 0
  diag(A) <- 0
  
  # Adding names
  colnames(A) <- rownames(A) <- gNames
  
  # Making the output a sparse matrix
  A <- as(A, 'dgCMatrix')
  
  # Return
  return(A)
}

makeNetworks <- function(X, nNet = 10, nCells = 500, nComp = 3, scaleScores = TRUE, symmetric = FALSE, q = 0.95){
  geneList <- rownames(X)
  nGenes <- length(geneList)
  nCol <- ncol(X)
  if(nGenes > 0){
    pbapply::pbsapply(seq_len(nNet), function(W){
      Z <- sample(x = seq_len(nCol), size = nCells, replace = TRUE)
      Z <- as.matrix(X[,Z])
      Z <- Z[apply(Z,1,sum) > 0,]
      if(nComp > 1 & nComp < nGenes){
        Z <- pcNet(Z, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = q, verbose = FALSE)  
      } else {
        stop('nComp should be greater or equal than 2 and lower than the total number of genes')
      }
      O <- matrix(data = 0, nrow = nGenes, ncol = nGenes)
      rownames(O) <- colnames(O) <- geneList
      O[rownames(Z), colnames(Z)] <- as.matrix(Z)
      O <- as(O, 'dgCMatrix')
      return(O)
    })  
  } else {
    stop('Gene names are required')
  }
}

## make networks
scTeni_makeNet <- function(X, norm_method = "new", nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                           nc_q = 0.05, td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5){
  
  if (norm_method == "new"){
    # Comparing gene ids.
    xNames <- rownames(X)
    nGenes <- length(xNames)
    
    # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
    set.seed(1)
    xList <- makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
  } else{
    # Comparing gene ids.
    xNames <- rownames(X)
    nGenes <- length(xNames)
    
    # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
    set.seed(1)
    xList <- scTenifoldNet::makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
  }
  oX <- list2network(xList)
  rownames(oX) <- xNames
  colnames(oX) <- xNames
  if (td_K == 0){
    return(oX)
  } else{
    set.seed(1)
    tensorOut <- scTenifoldNet::tensorDecomposition(xList, K = td_K, maxIter = td_maxIter, maxError = td_maxError)
    tx <- as.matrix(tensorOut$X)
    # Return
    return(tx)
  }
}
