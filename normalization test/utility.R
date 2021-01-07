##############################################
########## evaluation methods ################
##############################################

Acurracy_fun <- function(X, Q){
  X <- as.matrix(X)
  X <- X/max(abs(X))
  X <- (X - 0)
  X[abs(X) < quantile(abs(X), Q)] <-  0
  diag(X) <- 1
  dataMatrix <- X
  TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
  FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
  TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
  FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
  ACC <- round((TP+TN)/(TP+TN+FN+FP),2)
  return(ACC)
}

Recall_fun <- function(X, Q){
  X <- as.matrix(X)
  X <- X/max(abs(X))
  X <- (X - 0)
  X[abs(X) < quantile(abs(X), Q)] <-  0
  diag(X) <- 1
  dataMatrix <- X
  TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
  FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
  TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
  FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
  REC <- round((TP)/((40*40)+(58*58)),2)
  return(REC)
}

Fmeasure_fun <- function(X, Q){
  X <- as.matrix(X)
  X <- X/max(abs(X))
  X <- (X - 0)
  X[abs(X) < quantile(abs(X), Q)] <-  0
  diag(X) <- 1
  dataMatrix <- X
  TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
  FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
  TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
  FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
  Fmeasure <- round(TP / (TP + (FP + FN)/2),2)
  return(Fmeasure)
}

metricOutput_fun <- function(Acurracy, ReCall, Fmeasure, Q = 0, nCells){
  accMean <- apply(Acurracy, 2, mean)
  accSD <- apply(Acurracy, 2, sd)
  recallMean <- apply(ReCall, 2, mean)
  recallSD <- apply(ReCall, 2, sd)
  FmeasureMean <- apply(Fmeasure, 2, mean)
  FmeasureSD <- apply(Fmeasure, 2, sd)
  
  outputMetric <- NULL
  outputMetric$q <- rep(Q, length(nCells))
  outputMetric$nCells <- nCells
  outputMetric$accLB <- accMean-accSD
  outputMetric$acc <- accMean
  outputMetric$accUB <- accMean+accSD
  outputMetric$recallLB <- recallMean - recallSD
  outputMetric$recall <- recallMean
  outputMetric$recallUB <- recallMean + recallSD
  outputMetric$FmeasureLB <- FmeasureMean - FmeasureSD
  outputMetric$Fmeasure <- FmeasureMean
  outputMetric$FmeasureUB <- FmeasureMean + FmeasureSD
  
  outputMetric <- as.data.frame(outputMetric)
  outputMetric
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
