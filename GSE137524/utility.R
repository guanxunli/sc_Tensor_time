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


as.tensor <- function(x,drop=FALSE){
  stopifnot(is.array(x)||is.vector(x))
  tnsr <- list()
  if (is.vector(x)){
    modes <- c(length(x))
    num_modes <- 1L
    tnsr$modes <- modes
    tnsr$num_modes <- num_modes
    tnsr$data <- x
    #new("Tensor", num_modes, modes, data = x)
  }
  else {
    modes <- dim(x)
    num_modes <- length(modes)
    dim1s <- which(modes==1)
    if (drop && (length(dim1s)>0)){
      modes <- modes[-dim1s]
      num_modes <- num_modes-length(dim1s)
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- array(x,dim=modes)
      #new("list",num_modes,modes,data=array(x,dim=modes))
    } else {
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- x
    }
  }
  return(tnsr)
}

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

cpDecomposition <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5){
  kronecker_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for(i in 2:length(L)){
      retmat <- kronecker(retmat,L[[i]])
    }
    retmat
  }
  
  
  fnorm <- function(tnsr){
    arr<-tnsr$data
    sqrt(sum(arr*arr))
  }
  
  rs_unfold <- function(tnsr,m=NULL){
    if(is.null(m)) stop("mode m must be specified")
    num_modes <- tnsr$num_modes
    rs <- m
    cs <- (1:num_modes)[-m]
    unfold(tnsr,row_idx=rs,col_idx=cs)
  }
  
  unfold <- function(tnsr,row_idx=NULL,col_idx=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row and column indices must be specified")
    num_modes <- tnsr$num_modes
    if (length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
    if(any(rs<1) || any(rs>num_modes) || any(cs < 1) || any(cs>num_modes)) stop("illegal indices specified")
    perm <- c(rs,cs)
    if (any(sort(perm,decreasing=TRUE) != num_modes:1)) stop("missing and/or repeated indices")
    modes <- tnsr$modes
    mat <- tnsr$data
    new_modes <- c(prod(modes[rs]),prod(modes[cs]))
    #rearranges into a matrix
    mat <- aperm(mat,perm)
    dim(mat) <- new_modes
    as.tensor(mat)
  }
  
  hadamard_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for (i in 2:length(L)){
      retmat <- retmat*L[[i]]
    }
    retmat
  }
  
  khatri_rao_list <- function(L,reverse=FALSE){
    stopifnot(all(unlist(lapply(L,is.matrix))))
    ncols <- unlist(lapply(L,ncol))
    stopifnot(length(unique(ncols))==1)
    ncols <- ncols[1]
    nrows <- unlist(lapply(L,nrow))
    retmat <- matrix(0,nrow=prod(nrows),ncol=ncols)
    if (reverse) L <- rev(L)
    for(j in 1:ncols){
      Lj <- lapply(L,function(x) x[,j])
      retmat[,j] <- kronecker_list(Lj)
    }
    retmat
  }
  
  superdiagonal_tensor <- function(num_modes,len,elements=1L){
    modes <- rep(len,num_modes)
    arr <- array(0, dim = modes)
    if(length(elements)==1) elements <- rep(elements,len)
    for (i in 1:len){
      txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
      eval(parse(text=txt))
    }
    as.tensor(arr)
  }
  
  ttl<-function(tnsr,list_mat,ms=NULL){
    if(is.null(ms)||!is.vector(ms)) stop ("m modes must be specified as a vector")
    if(length(ms)!=length(list_mat)) stop("m modes length does not match list_mat length")
    num_mats <- length(list_mat)
    if(length(unique(ms))!=num_mats) warning("consider pre-multiplying matrices for the same m for speed")
    mat_nrows <- vector("list", num_mats)
    mat_ncols <- vector("list", num_mats)
    for(i in 1:num_mats){
      mat <- list_mat[[i]]
      m <- ms[i]
      mat_dims <- dim(mat)
      modes_in <- tnsr$modes
      stopifnot(modes_in[m]==mat_dims[2])
      modes_out <- modes_in
      modes_out[m] <- mat_dims[1]
      tnsr_m <- rs_unfold(tnsr,m=m)$data
      retarr_m <- mat%*%tnsr_m
      tnsr <- rs_fold(retarr_m,m=m,modes=modes_out)
    }	
    tnsr	
  }
  
  rs_fold <- function(mat,m=NULL,modes=NULL){
    if(is.null(m)) stop("mode m must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    num_modes <- length(modes)
    rs <- m
    cs <- (1:num_modes)[-m]
    fold(mat,row_idx=rs,col_idx=cs,modes=modes)
  }
  
  fold <- function(mat, row_idx = NULL, col_idx = NULL, modes=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    if(!is(mat,"list")){
      if(!is.matrix(mat))  stop("mat must be of class 'matrix'")
    }else{
      stopifnot(mat$num_modes==2)
      mat <- mat$data			
    }
    num_modes <- length(modes)
    stopifnot(num_modes==length(rs)+length(cs))
    mat_modes <- dim(mat)
    if((mat_modes[1]!=prod(modes[rs])) || (mat_modes[2]!=prod(modes[cs]))) stop("matrix nrow/ncol does not match Tensor modes")
    #rearranges into array
    iperm <- match(1:num_modes,c(rs,cs))
    as.tensor(aperm(array(mat,dim=c(modes[rs],modes[cs])),iperm))
  }
  
  
  if(is.null(num_components)) stop("num_components must be specified")
  stopifnot(is(tnsr,"list"))
  #if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
  
  #initialization via truncated hosvd
  num_modes <- tnsr$num_modes
  modes <- tnsr$modes
  U_list <- vector("list",num_modes)
  unfolded_mat <- vector("list",num_modes)
  tnsr_norm <- fnorm(tnsr)
  for(m in 1:num_modes){
    unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)$data
    U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
  }
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  #set up convergence check
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(est){
    curr_resid <- fnorm(as.tensor(est$data - tnsr$data))
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter==1) return(FALSE)
    if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
    else{ return(FALSE)}
  }	
  #progress bar
  pb <- txtProgressBar(min=0,max=max_iter,style=3)
  #main loop (until convergence or max_iter)
  norm_vec <- function(vec){
    norm(as.matrix(vec))
  }
  while((curr_iter < max_iter) && (!converged)){
    setTxtProgressBar(pb,curr_iter)
    for(m in 1:num_modes){
      V <- hadamard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
      V_inv <- solve(V)			
      tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
      lambdas <- apply(tmp,2,norm_vec)
      U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
      Z <- superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
      est <- ttl(Z,U_list,ms=1:num_modes)
    }
    #checks convergence
    if(CHECK_CONV(est)){
      converged <- TRUE
      setTxtProgressBar(pb,max_iter)
    }else{
      curr_iter <- curr_iter + 1
    }
  }
  if(!converged){setTxtProgressBar(pb,max_iter)}
  close(pb)
  #end of main loop
  #put together return list, and returns
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]
  norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
  invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
}


# makeNetworks <- function(X, nNet = 10, nCells = 500, nComp = 3, scaleScores = TRUE, symmetric = FALSE, q = 0.95){
#   geneList <- rownames(X)
#   nGenes <- length(geneList)
#   nCol <- ncol(X)
#   if(nGenes > 0){
#     pbapply::pbsapply(seq_len(nNet), function(W){
#       Z <- sample(x = seq_len(nCol), size = nCells, replace = TRUE)
#       Z <- as.matrix(X[,Z])
#       Z <- Z[apply(Z,1,sum) > 0,]
#       if(nComp > 1 & nComp < nGenes){
#         Z <- pcNet(Z, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = q, verbose = FALSE)  
#       } else {
#         stop('nComp should be greater or equal than 2 and lower than the total number of genes')
#       }
#       O <- matrix(data = 0, nrow = nGenes, ncol = nGenes)
#       rownames(O) <- colnames(O) <- geneList
#       O[rownames(Z), colnames(Z)] <- as.matrix(Z)
#       O <- as(O, 'dgCMatrix')
#       return(O)
#     })  
#   } else {
#     stop('Gene names are required')
#   }
# }

# ## make networks
# scTeni_makeNet <- function(X, norm_method = "new", nc_nNet = 10, nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
#                            nc_q = 0.05, td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5){

#   if (norm_method == "new"){
#     # Comparing gene ids.
#     xNames <- rownames(X)
#     nGenes <- length(xNames)

#     # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
#     set.seed(1)
#     xList <- makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
#   } else{
#     # Comparing gene ids.
#     xNames <- rownames(X)
#     nGenes <- length(xNames)

#     # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
#     set.seed(1)
#     xList <- scTenifoldNet::makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
#   }
#   oX <- list2network(xList)
#   rownames(oX) <- xNames
#   colnames(oX) <- xNames
#   if (td_K == 0){
#     return(oX)
#   } else{
#     set.seed(1)
#     tensorOut <- scTenifoldNet::tensorDecomposition(xList, K = td_K, maxIter = td_maxIter, maxError = td_maxError)
#     tx <- as.matrix(tensorOut$X)
#     # Return
#     return(tx)
#   }
# }

######################################################
############ Other regression function ###############
######################################################
my_regression <- function(network_list, time_vec) {
  ## Getting beta mat
  n_net <- length(network_list)
  X_mat <- cbind(1, time_vec)
  tmp <- as.numeric(network_list[[1]])
  nGenes <- nrow(network_list[[1]])
  Y <- tmp
  for (i in 2:n_net){
    tmp <- as.numeric(network_list[[i]])
    Y <- rbind(Y, tmp)
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

UMAP_order <- function(dta){
  dta_svd <- RSpectra::svds(dta, k = 50)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  UMAP <- umap(dta_pcscore)
  rownames(UMAP) <- rownames(dta)
  ## GSEA
  mDistance <- mahalanobis(UMAP, colMeans(UMAP), cov(UMAP))
  BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
  set.seed(1)
  E <- fgseaMultilevel(BIOP, mDistance[!grepl('^RPL|^RPS|^RP[[:digit:]]|^MT-',names(mDistance))])
  E <- E[order(E$pval, decreasing = FALSE),]
  E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
  return(E)
}