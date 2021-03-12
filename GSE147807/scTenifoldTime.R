source("pcNet.R")
source("tensor_utility.R")
source("tensorDecomposition.R")
source("regression_fun.R")

# dta_list is the list for dta
# time_vec is the time vector for data list
# nComp is # of pc for PC regression
# q is the quntitle for pcNet
# K is # of PC in tensor decomposition
scTenifoldTime_beta <- function(dta_list, method = "pcnet", time_vec, nComp = 5, q = 0,
                           K = 10, maxIter = 10000, maxError = 1e-5){
  res <- list()
  n_net <- length(dta_list)
  nGenes <- nrow(dta_list[[1]])
  gene_name <- rownames(dta_list[[1]])
  ## generate network
  network_list <- list()
  for (i in seq_len(n_net)) {
    if (method == "pcnet"){
      network_list[[i]] <- pcNet(dta_list[[i]], nComp = nComp, scaleScores = TRUE, symmetric = FALSE, q = q, verbose = TRUE)
      network_list[[i]] <- round(network_list[[i]], 2)
    } else{
      network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
      network_list[[i]] <- round(network_list[[i]], 2)
    }
    rownames(network_list[[i]]) <- gene_name
    print(paste0("Finish network ", i))
  }
  res$network_list <- network_list
  ## tensor decomposition
  network_tensor <- tensorDecomposition(network_list, K = K, maxIter = maxIter, maxError = maxError)
  for (i in seq_len(length(network_tensor))) {
    diag(network_tensor[[i]]) <- 1
  }
  res$network_tensor <- network_tensor
  print("Finish tensor decomposition part.")

  ## Do regression
  res_regression <- my_regression(network_list = network_list, time_vec = time_vec)
  beta_mat <- res_regression$beta_mat
  rownames(beta_mat) <- colnames(beta_mat) <- rownames(dta_list[[1]])
  t_mat <- res_regression$t_mat
  rm(res_regression)
  res$beta_mat <- beta_mat
  res$t_mat <- t_mat
  
  ## Do regression
  res_regression <- my_regression(network_list = network_tensor, time_vec = time_vec)
  beta_mat <- res_regression$beta_mat
  rownames(beta_mat) <- colnames(beta_mat) <- rownames(dta_list[[1]])
  t_mat <- res_regression$t_mat
  rm(res_regression)
  res$beta_mat_tensor <- beta_mat
  res$t_mat_tensor <- t_mat

  ## return results
  return(res)
}

scTenifoldTime_tensor <- function(dta_list, time_vec, method = "pcnet", nComp = 5, q = 0,
                           K = 10, maxIter = 10000, maxError = 1e-5, thres = 0.05, nDecimal = 2,
                           ma_nDim = 3, scoreType = "pos", eps = 0){
  res <- list()
  n_net <- length(dta_list)
  nGenes <- nrow(dta_list[[1]])
  gene_name <- rownames(dta_list[[1]])
  ## generate network
  network_list <- list()
  for (i in seq_len(n_net)) {
    if (method == "pcnet"){
      network_list[[i]] <- pcNet(dta_list[[i]], nComp = nComp, scaleScores = TRUE, symmetric = FALSE, q = q, verbose = TRUE)
      network_list[[i]] <- round(network_list[[i]], 2)
    } else{
      network_list[[i]] <- cor(t(dta_list[[i]]), method = "spearman")
      network_list[[i]] <- round(network_list[[i]], 2)
    }
    rownames(network_list[[i]]) <- gene_name
    print(paste0("Finish network ", i))
  }
  res$network_list <- network_list
  ## tensor decomposition
  set.seed(1)
  tensor_output <- tensorDecomposition_time(network_list, K = K, maxError = maxError, maxIter = maxIter, 
                                            time_vec = time_vec, thres = thres, nDecimal = nDecimal)
  res$tensor_output <- tensor_output
  print("Finish tensor decomposition part.")
  
  if (is.null(tensor_output$network0) == TRUE){
    return(res)
  } else{
    ## Split of tensor output
    tX <- as.matrix(tensor_output$network0)
    tY <- as.matrix(tensor_output$network1)
    rownames(tX) <- colnames(tX) <- gene_name
    rownames(tY) <- colnames(tY) <- gene_name
    
    ## Making it symmetric to fulfill the requirements of the MA
    tX <- (tX + t(tX))/2
    tY <- (tY + t(tY))/2
    
    ## Non-linear manifold alignment
    set.seed(1)
    mA <- scTenifoldNet::manifoldAlignment(tX , tY, d = ma_nDim)
    res$mA <- mA
    mA_X <- mA[1:nGenes, ]
    mA_Y <- mA[-(1:nGenes), ]
    print("Finish manifold alignment part.")
    
    ## return order of gene
    gene_diff <- rowSums((mA_X - mA_Y)^2)
    names(gene_diff) <- gene_name
    gene_diff <- sort(gene_diff, decreasing = TRUE)
    res$gene_diff <- gene_diff
    
    return(res)
  }
}