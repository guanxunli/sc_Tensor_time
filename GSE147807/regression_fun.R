my_regression <- function(network_list, time_vec) {
  ## Getting beta mat
  n_net <- length(network_list)
  X_mat <- cbind(1, time_vec)
  tmp <- as.numeric(network_list[[1]])
  nGenes <- nrow(network_list[[1]])
  Y <- tmp
  if (n_net > 1){
    for (i in 2:n_net){
      tmp <- as.numeric(network_list[[i]])
      Y <- rbind(Y, tmp)
    }
  }
  ## Do regress
  beta_coef <- solve(crossprod(X_mat), crossprod(X_mat, Y))
  beta_mat <- round(matrix(beta_coef[2, ], nrow = nGenes), 2)
  diag(beta_mat) <- 1
  beta_adj <- abs(beta_mat)
  beta_adj[which(beta_adj > 0)] <- 1
  
  ## t-test
  lm_t_p <- matrix(1, nrow = ncol(X_mat) - 1, ncol = dim(Y)[2])
  df <- nrow(X_mat) - ncol(X_mat)
  index <- setdiff(which(beta_coef[1, ] != 0), which(beta_coef[2, ] == 0))
  Y_fit <- X_mat %*% beta_coef[, index]
  sigma_fit <- colSums((Y[, index] - Y_fit)^2) / df
  
  for (i in 2:ncol(X_mat)){
    sd_t <- sqrt(solve(crossprod(X_mat))[i, i] * sigma_fit)
    t_test <- beta_coef[i, index] / sd_t
    t_test_p <- 2 * (1 - pt(abs(t_test), df = df))
    lm_t_p[i - 1, index] <- t_test_p
  }
  t_mat <- matrix(lm_t_p, nrow = nGenes)
  
  ## return results
  res <- list()
  res$beta_mat <- beta_mat
  res$t_mat <- t_mat
  res$beta_adj <- beta_adj
  return(res)
}