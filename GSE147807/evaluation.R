gene_list_fun <- function(res, ncomp = 50) {
  gene_list <- list()
  
  #### Without tensor decomposition
  beta_mat <- res$beta_mat
  t_mat <- res$t_mat
  ## with t-test filter
  beta_use <- matrix(0, nrow = nrow(beta_mat), ncol = ncol(beta_mat)) 
  beta_use[which(t_mat < 0.05)] <- beta_mat[which(t_mat < 0.05)] 
  rownames(beta_use) <- rownames(beta_mat) 
  colnames(beta_use) <- colnames(beta_mat)
  
  ## PCA without t-test
  set.seed(1)
  dta_svd <- RSpectra::svds(beta_mat, k = ncomp)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  # UMAP
  UMAP <- uwot::umap(dta_pcscore)
  rownames(UMAP) <- rownames(beta_mat)
  mDistance <- mahalanobis(UMAP, colMeans(UMAP), cov(UMAP))
  names(mDistance) <- toupper(rownames(beta_mat))
  mDistance <- sort(mDistance, decreasing = TRUE)
  gene_list$UMAP <- mDistance
  # TSNE
  TSNE <- Rtsne::Rtsne(dta_pcscore, pca = FALSE, check_duplicates = FALSE)
  D <- mahalanobis(TSNE$Y, center = colMeans(TSNE$Y), cov = cov(TSNE$Y))
  names(D) <- toupper(rownames(beta_mat))
  D <- sort(D, decreasing = TRUE)
  gene_list$TSNE <- D
  
  ## PCA with t-test filter
  set.seed(1)
  dta_svd <- RSpectra::svds(beta_use, k = ncomp)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  # UMAP
  UMAP <- uwot::umap(dta_pcscore)
  rownames(UMAP) <- rownames(beta_use)
  mDistance <- mahalanobis(UMAP, colMeans(UMAP), cov(UMAP))
  names(mDistance) <- toupper(rownames(beta_use))
  mDistance <- sort(mDistance, decreasing = TRUE)
  gene_list$UMAP_filter <- mDistance
  # TSNE
  TSNE <- Rtsne::Rtsne(dta_pcscore, pca = FALSE, check_duplicates = FALSE)
  D <- mahalanobis(TSNE$Y, center = colMeans(TSNE$Y), cov = cov(TSNE$Y))
  names(D) <- toupper(rownames(beta_use))
  D <- sort(D, decreasing = TRUE)
  gene_list$TSNE_filter <- D
  
  #### With tensor decomposition
  beta_mat <- res$beta_mat_tensor
  t_mat <- res$t_mat_tensor
  ## with t-test filter
  beta_use <- matrix(0, nrow = nrow(beta_mat), ncol = ncol(beta_mat)) 
  beta_use[which(t_mat < 0.05)] <- beta_mat[which(t_mat < 0.05)] 
  rownames(beta_use) <- rownames(beta_mat) 
  colnames(beta_use) <- colnames(beta_mat)
  
  ## PCA without t-test
  set.seed(1)
  dta_svd <- RSpectra::svds(beta_mat, k = ncomp)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  # UMAP
  UMAP <- uwot::umap(dta_pcscore)
  rownames(UMAP) <- rownames(beta_mat)
  mDistance <- mahalanobis(UMAP, colMeans(UMAP), cov(UMAP))
  names(mDistance) <- toupper(rownames(beta_mat))
  mDistance <- sort(mDistance, decreasing = TRUE)
  gene_list$UMAP_tensor <- mDistance
  # TSNE
  TSNE <- Rtsne::Rtsne(dta_pcscore, pca = FALSE, check_duplicates = FALSE)
  D <- mahalanobis(TSNE$Y, center = colMeans(TSNE$Y), cov = cov(TSNE$Y))
  names(D) <- toupper(rownames(beta_mat))
  D <- sort(D, decreasing = TRUE)
  gene_list$TSNE_tensor <- D
  
  ## PCA with t-test filter
  set.seed(1)
  dta_svd <- RSpectra::svds(beta_use, k = ncomp)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  # UMAP
  UMAP <- uwot::umap(dta_pcscore)
  rownames(UMAP) <- rownames(beta_use)
  mDistance <- mahalanobis(UMAP, colMeans(UMAP), cov(UMAP))
  names(mDistance) <- toupper(rownames(beta_use))
  mDistance <- sort(mDistance, decreasing = TRUE)
  gene_list$UMAP_tensor_filter <- mDistance
  # TSNE
  TSNE <- Rtsne::Rtsne(dta_pcscore, pca = FALSE, check_duplicates = FALSE)
  D <- mahalanobis(TSNE$Y, center = colMeans(TSNE$Y), cov = cov(TSNE$Y))
  names(D) <- toupper(rownames(beta_use))
  D <- sort(D, decreasing = TRUE)
  gene_list$TSNE_tensor_filter <- D
  
  ## return results
  return(gene_list)
}

gene_rep_fun <- function(res) {
  index1 <- res$tensor_output$index1
  n_index <- length(index1)
  gene_rep <- matrix(NA, nrow = nrow(res$network_list[[1]]), ncol = 2 * n_index)
  rownames(gene_rep) <- rownames(res$network_list[[1]])
  lambda <- matrix(res$tensor_output$lambdas[index1], nrow = nrow(gene_rep), 
                   ncol = n_index, byrow = TRUE)
  gene_rep[ ,seq_len(n_index)] <- res$tensor_output$U[[1]][, index1] * lambda
  gene_rep[ ,(n_index + 1):(2 * n_index)] <- res$tensor_output$U[[2]][, index1] * lambda
  return(gene_rep)
}

gene_list_fun_tensor <- function(res) {
  gene_list <- list()
  
  ## Compare two networks
  gene_diff <- res$tensor_output$gene_diff
  gene_diff <- sort(gene_diff, decreasing = TRUE)
  gene_list$network <- gene_diff
  
  ## After manifold alignment
  nGenes <- nrow(res$network_list[[1]])
  gene_name <- rownames(res$network_list[[1]])
  iter <- 2
  for (i in c(2,3,5,10,20,30)) {
    mA_tmp <- res$mA[, 1:i]
    mA_X <- mA_tmp[1:nGenes, ]
    mA_Y <- mA_tmp[-(1:nGenes), ]
    ## return order of gene
    gene_diff <- rowSums(abs((mA_X - mA_Y)))
    names(gene_diff) <- gene_name
    gene_diff <- sort(gene_diff, decreasing = TRUE)
    gene_list[[iter]] <- gene_diff
    names(gene_list)[iter] <- paste0("manifoldalignment", i)
    iter <- iter + 1
  }
  
  ## Tensor representation
  gene_rep <- gene_rep_fun(res)
  mDistance <- mahalanobis(gene_rep, colMeans(gene_rep), cov(gene_rep))
  names(mDistance) <- toupper(rownames(gene_rep))
  mDistance <- sort(mDistance, decreasing = TRUE)
  gene_list$tensor_rep <- mDistance
  
  ## return results
  return(gene_list)
}