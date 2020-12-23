library(Matrix)
library(ComplexHeatmap)
library(ggplot2)

dta <- read.csv("SERGIO/diffOutput_S_1000.csv")
range(dta)
dim(dta)
rownames(dta) <- paste0("Gene", 1:nrow(dta))
colnames(dta) <- paste0("Cell", 1:ncol(dta))
source("code/utility.R")
nGenes <- nrow(dta)
#############################################################################
############## Using X_{t} to predict X_{t} #################################
#############################################################################

# png(paste0('network_',nComp,'q_', q, 'nCells_', nCells, '.png'), width = 1000, height = 1000, res = 300)
# print(Heatmap(dta_net, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
# dev.off()

#### All data sets
## old method
dta_net <- scTeni_makeNet(dta, norm_method = "old", nc_q = 0.3, tensor_dc = "FALSE", nc_nCells = 800)
Heatmap(dta_net, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
dta_net <- scTeni_makeNet(dta, norm_method = "old", nc_q = 0.3, tensor_dc = "TRUE", nc_nCells = 800, td_K = 20)
Heatmap(dta_net, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
## new method
dta_net <- scTeni_makeNet(dta, norm_method = "new", nc_q = 0.3, tensor_dc = "FALSE", nc_nCells = 800)
Heatmap(dta_net, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
dta_net <- scTeni_makeNet(dta, norm_method = "new", nc_q = 0.3, tensor_dc = "TRUE", nc_nCells = 800, td_K = 20)
Heatmap(dta_net, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)

#### 10 data set
## divide into different data sets
scRNAsq_list <- list()
for (i in 1:10){
  scRNAsq_list[[i]] <- dta[, (1 + 1000 * (i - 1)):(1000 * i)]
}

## generate old networks
Network_list_old <- list()
for (i in 1:10){
  Network_list_old[[i]] <- scTeni_makeNet(scRNAsq_list[[i]], norm_method = "old", nc_q = 0.45, tensor_dc = "FALSE", nc_nCells = 800, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('oMatrix_old',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_old[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_old)
saveRDS(reg_res, "reg_res_old.rds")

Network_list_old_td <- list()
for (i in 1:10){
  Network_list_old_td[[i]] <- scTeni_makeNet(scRNAsq_list[[i]],norm_method = "old", nc_q = 0.45, tensor_dc = "TRUE", nc_nCells = 800, td_K = 20, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('tMatrix_old',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_old_td[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_old_td)
saveRDS(reg_res, "reg_res_old_td.rds")

## generate new networks
Network_list_new <- list()
for (i in 1:10){
  Network_list_new[[i]] <- scTeni_makeNet(scRNAsq_list[[i]], norm_method = "new", nc_q = 0.45, tensor_dc = "FALSE", nc_nCells = 800, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('oMatrix_new',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_new[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_new)
saveRDS(reg_res, "reg_res_new.rds")

Network_list_new_td <- list()
for (i in 1:10){
  Network_list_new_td[[i]] <- scTeni_makeNet(scRNAsq_list[[i]],norm_method = "new", nc_q = 0.45, tensor_dc = "TRUE", nc_nCells = 800, td_K = 20, nc_nComp = 5)
}

for(i in 1:10){
  png(paste0('tMatrix_new',i,'.png'), width = 1000, height = 1000, res = 300)
  print(Heatmap(Network_list_new_td[[i]], row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  dev.off()
}

reg_res <- regression_res(Network_list_new_td)
saveRDS(reg_res, "reg_res_new_td.rds")

#### Community detection
plot_coef <- function(x_rds){
  reg_beta <- x_rds$coef
  reg_beta <- matrix(reg_beta[2, ], nrow = nGenes, ncol = nGenes)
  print(Heatmap(reg_beta, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE))
  reg_beta_pca <- prcomp(reg_beta, scale. = FALSE)
  print(ggplot(data = NULL, aes(x = reg_beta_pca$x[, 1], y = reg_beta_pca$x[, 2])) +
    geom_point(color = as.factor(c(rep(1, 28), rep(2, 21), rep(3, 30), rep(4, 20), rep(5, 3)))) +
    labs(x = "PC1", y = "PC2", title = "PCA results"))
  reg_beta_tsne <- Rtsne::Rtsne(reg_beta)
  reg_beta_tsne <- reg_beta_tsne$Y
  print(ggplot(data = NULL, aes(x = reg_beta_tsne[, 1], y = reg_beta_tsne[, 2])) +
          geom_point(color = as.factor(c(rep(1, 28), rep(2, 21), rep(3, 30), rep(4, 20), rep(5, 3)))) +
          labs(x = "tsne1", y = "tsne2", title = "tsne results"))
}
reg_res <- readRDS("results/simulation/reg_res_old.rds")
plot_coef(reg_res)
reg_res <- readRDS("results/simulation/reg_res_old_td.rds")
plot_coef(reg_res)
reg_res <- readRDS("results/simulation/reg_res_new.rds")
plot_coef(reg_res)
reg_res <- readRDS("results/simulation/reg_res_new_td.rds")
plot_coef(reg_res)


#############################################################################
############## Using X_{t-1} to predict X_{t} ###############################
#############################################################################
# divide into different data sets
scRNAsq_list <- list()
for (i in 1:10){
  scRNAsq_list[[i]] <- dta[, (1 + 1000 * (i - 1)):(1000 * i)]
}

pcNet <- function(X, Y, nComp = 3, scaleScores = TRUE,symmetric = FALSE, q = 0, verbose = TRUE) {
  if (!all(Matrix::rowSums(X) > 0)) {
    stop('Quality control has not been applied over the matrix.')
  }
  xClass <- class(X)[[1]]
  validClass <- xClass %in% c('matrix', 'dgCMatrix')
  if (!validClass) {
    stop('Input should be a matrix with cells as columns and genes as rows')
  }
  if (nComp < 2 | nComp >= nrow(X)) {
    stop('nCom should be greater or equal than 2 and lower than the total number of genes')
  }
  if (!identical(dim(X), dim(Y))){
    stop('Dimension of two matrices are not equal.')
  }
  gNames <- rownames(X)
  pcCoefficients <- function(K) {
    # Taking out the gene to be regressed out
    Xi <- X
    Xi[, K] <- Y[, K]
    Xi <- scale(Xi)
    y <- Xi[, K]
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
  
  # Standardizing the data
  # X <- (scale(Matrix::t(X)))
  X <- Matrix::t(X)
  Y <- Matrix::t(Y)
  
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

makeNetworks <- function(X, Y, nNet = 10, nCells = 500, nComp = 3, scaleScores = TRUE, symmetric = FALSE, q = 0.95){
  geneList <- rownames(X)
  nGenes <- length(geneList)
  nCol <- ncol(X)
  if(nGenes > 0){
    pbapply::pbsapply(seq_len(nNet), function(W){
      cell_index <- sample(x = seq_len(nCol), size = nCells, replace = TRUE)
      ZX <- as.matrix(X[, cell_index])
      ZY <- as.matrix(Y[, cell_index])
      intersect_index <- intersect(which(rowSums(ZX) > 0), which(rowSums(ZY) > 0))
      ZX <- ZX[intersect_index, ]
      ZY <- ZY[intersect_index, ]
      if(nComp > 1 & nComp < nGenes){
        Z <- pcNet(ZX, ZY, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = q, verbose = FALSE)
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

generate_network <- function(dta_list, nComp = 3, q = 0.95, nCells = 500){
  n_list <- length(dta_list)
  for (iter in 1:(n_list - 1)){
    dta_net_list <- makeNetworks(X = dta_list[[iter]], Y = dta_list[[iter + 1]], nComp = nComp, q = q, nCells = nCells)
    dta_net <- scTenifoldNet::tensorDecomposition(xList = dta_net_list)
    dta_net <- as.matrix(dta_net$X)
    rownames(dta_net) <- NULL
    colnames(dta_net) <- NULL
    nGenes <- nrow(dta)
    png(paste0('network_',iter, 'nComp_', nComp, 'q_', q, 'nCells_', nCells, '.png'), width = 1000, height = 1000, res = 300)
    Heatmap(dta_net, row_order = seq_len(nGenes), column_order = seq_len(nGenes), show_heatmap_legend = FALSE)
    dev.off()
  }
}

generate_network(scRNAsq_list, nComp = 5, q = 0.85, nCells = 800)

