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