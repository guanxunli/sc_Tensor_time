source("tensor_utility.R")

tensorDecomposition <- function(xList, K = 5, maxError = 1e-5, maxIter = 1e3){
  xNets <- length(xList)
  nNet <- xNets
  xGenes <- unique(unlist(lapply(xList, rownames)))
  sGenes <- xGenes
  nGenes <- length(sGenes)

  tensorX <- array(data = 0, dim = c(nGenes,nGenes,1,nNet))

  for(i in seq_len(nNet)){
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
    tensorX[,,,i] <- tempX
  }

  set.seed(1)
  tensorX <- as.tensor(tensorX)
  tensorX <- cpDecomposition(tnsr = tensorX, num_components = K, max_iter = maxIter, tol = maxError)
  tensorOutput <- list()
  for (i in seq_len(nNet)){
    tensorOutput[[i]] <- round(tensorX$est$data[,,,i], 2)
  }
  return(tensorOutput)
}

# time_vec is the time vector for network list
# thres is the thresdhold for tensordecomposition
tensorDecomposition_time <- function(xList, K = 5, maxError = 1e-5, maxIter = 1e3, time_vec, thres = 0.05, nDecimal = 2){
  xNets <- length(xList)
  nNet <- xNets
  xGenes <- unique(unlist(lapply(xList, rownames)))
  sGenes <- xGenes
  nGenes <- length(sGenes) 
  
  tensorX <- array(data = 0, dim = c(nGenes,nGenes, nNet))
  
  for(i in seq_len(nNet)){
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
    tensorX[,,i] <- tempX
  }
  
  set.seed(1)
  tensorX <- as.tensor(tensorX)
  tensorX <- cpDecomposition(tnsr = tensorX, num_components = K, max_iter = maxIter, tol = maxError)
  tensor_time <- tensorX$U[[3]]
  
  ## Do regression
  X_mat <- cbind(1, time_vec)
  Y <- tensor_time
  beta_coef <- solve(crossprod(X_mat), crossprod(X_mat, Y))
    
  ## t-test
  lm_t_p <- rep(NA, dim(Y)[2])
  df <- nrow(X_mat) - ncol(X_mat)
  Y_fit <- X_mat %*% beta_coef
  sigma_fit <- colSums((Y - Y_fit)^2) / df
  sd_t <- sqrt(solve(crossprod(X_mat))[2, 2] * sigma_fit)
  t_test <- beta_coef[2, ] / sd_t
  t_test_p <- 2 * (1 - pt(abs(t_test), df = df))
  
  ## Based on the t-value, get two tensor
  index1 <- which(t_test_p < thres)
  index0 <- setdiff(seq_len(length(t_test_p)), index1)
  
  if (length(index0) == 0){
    tensorOutput <- list()
    tensorOutput$lambdas <- tensorX$lambdas
    tensorOutput$U <- tensorX$U
    return(tensorOutput)
    print("All components are varied during the time, please set a smaller threshold for p-value.")
  } else if (length(index1) == 0) {
    tensorOutput <- list()
    tensorOutput$lambdas <- tensorX$lambdas
    tensorOutput$U <- tensorX$U
    return(tensorOutput)
    print('All components are not varied during the time, please set a bigger threshold for p-value.')
  } else{
    tensor0 <- array(data = 0, dim = c(nGenes,nGenes, nNet))
    for (i in seq_len(length(index0))) {
      index_use <- index0[i]
      tmp <- tensorX$lambdas[index_use] * tcrossprod(tensorX$U[[1]][, index_use], tensorX$U[[2]][, index_use])
      for (j in seq_len(nNet)){
        tensor0[,,j] <- tensor0[,,j] + tmp * tensorX$U[[3]][j, index_use]
      }
    }
    
    tensor1 <- array(data = 0, dim = c(nGenes,nGenes, nNet))
    for (i in seq_len(length(index1))) {
      index_use <- index1[i]
      tmp <- tensorX$lambdas[index_use] * tcrossprod(tensorX$U[[1]][, index_use], tensorX$U[[2]][, index_use])
      for (j in seq_len(nNet)){
        tensor1[,,j] <- tensor1[,,j] + tmp * tensorX$U[[3]][j, index_use]
      }
    }
  }
  
  ## Transfer the tensor to the network
  network0 <- matrix(0, nrow = nGenes, ncol = nGenes)
  network1 <- matrix(0, nrow = nGenes, ncol = nGenes)
  for (i in seq_len(nNet)) {
    network0 <- network0 + tensor0[,,i]
    network1 <- network1 + tensor1[,,i]
  }
  network0 <- network0 / nNet
  network1 <- network1 / nNet
  
  network0 <- round(network0, nDecimal)
  network1 <- round(network1, nDecimal)
  
  network0 <- as(network0, 'dgCMatrix')
  network1 <- as(network1, 'dgCMatrix')
  
  rownames(network0) <- colnames(network0) <- sGenes
  rownames(network1) <- colnames(network1) <- sGenes
  
  ## return results
  tensorOutput <- list()
  tensorOutput$lambdas <- tensorX$lambdas
  tensorOutput$U <- tensorX$U
  tensorOutput$network0 <- network0
  tensorOutput$network1 <- network1
  return(tensorOutput)
}