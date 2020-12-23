library(RSpectra)
library(apcluster)
library(igraph)
library(FNN)
library(NMF)

# generate test matrix
test_mat <- matrix(0, nrow = 10, ncol = 10)
test_mat[1:5, 1:5] <- matrix(abs(runif(25)), nrow = 5)
test_mat[6:10, 6:10] <- matrix(abs(runif(25)), nrow = 5)
diag(test_mat[1:5, 6:10]) <- runif(5, max = 0.1)
diag(test_mat[6:10, 1:5]) <- runif(5, max = 0.1)
test_mat

#########################################
#### Affinity Propagation Clustering ####
#########################################
APC_fun <- function(net){
  res_ap <- apcluster(s = negDistMat(), x = net)
  com_det_ap <- res_ap@clusters
  res <- rep(NA, nrow(net))
  for (i in 1:length(com_det_ap)){
    res[as.numeric(com_det_ap[[i]])] <- i
  }
  return(res)
}

APC_res <- APC_fun(test_mat)

#########################################
#### SCORE method  ######################
#########################################
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
  SCORE_kmeans <- kmeans(eig_ratio, centers = K, nstart = 10)
  return(SCORE_kmeans$cluster)
}

res_SCORE <- SCORE_fun(net = test_mat, K = 2)

#########################################
#### label propagation algorithm (LPA) ##
#########################################
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

LPA_res <- LPA_fun(test_mat, K = 3)

################################################
#### non-negative matrix factorlization (NMF) ##
################################################
# net is the network
# K is the number of communities
NMF_fun <- function(net, K){
  net <- abs(net)
  res_NMF <- nmf(x = net, rank = K)@fit@W
  res_NMF <- apply(res_NMF, 1, which.max)
  return(res_NMF)
}

NMF_res <- NMF_fun(test_mat, K = 2)


