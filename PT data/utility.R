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
  
  X <- scTenifoldNet::scQC(X, minLibSize = 1000, removeOutlierCells = TRUE, minPCT = 0.05, maxMTratio = 0.1)
  
  # Counts per million (CPM) normalization
  if (norm_method == "new"){
    X <- new_Normalization(X)
    # Comparing gene ids.
    xNames <- rownames(X)
    nGenes <- length(xNames)
    
    # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
    set.seed(1)
    xList <- makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
  } else{
    X <- scTenifoldNet::cpmNormalization(X)
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
  set.seed(1)
  tensorOut <- scTenifoldNet::tensorDecomposition(xList, K = td_K, maxIter = td_maxIter, maxError = td_maxError)
  res <- list()
  res$oX <- oX
  res$td <- tensorOut
  # Return
  return(res)
}