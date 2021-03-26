library(fgsea)
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')

## Given gene order
pathway_order <- function(gene_diff, scoreType = "pos", eps = 0) {
  set.seed(1)
  names(gene_diff) <- toupper(names(gene_diff))
  E <- fgseaMultilevel(BIOP, gene_diff[!grepl('^RPL|^RPS|^RP[[:digit:]]|^MT-',names(gene_diff))],
                       scoreType = scoreType, eps = eps)
  E <- E[order(E$pval, decreasing = FALSE),]
  E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
  return(E)
}

## Given beta matrix 
UMAP_order <- function(dta, ncomp = 50, scoreType = "pos", eps = 0) {
  ## Do UMAP
  set.seed(1)
  dta_svd <- RSpectra::svds(dta, k = ncomp)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  UMAP <- uwot::umap(dta_pcscore)
  rownames(UMAP) <- rownames(dta)
  ## GSEA
  mDistance <- mahalanobis(UMAP, colMeans(UMAP), cov(UMAP))
  names(mDistance) <- toupper(rownames(dta))
  set.seed(1)
  E <- pathway_order(gene_diff = mDistance, scoreType = scoreType, eps = eps)
  return(E)
}

TSNE_order <- function(dta, ncomp = 50, scoreType = "pos", eps = 0) {
  ## Do TSNE
  set.seed(1)
  dta_svd <- RSpectra::svds(dta, k = ncomp)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  TSNE <- Rtsne::Rtsne(dta_pcscore)
  ## GSEA
  D <- mahalanobis(TSNE$Y, center = colMeans(TSNE$Y), cov = cov(TSNE$Y))
  names(D) <- toupper(rownames(dta))
  set.seed(1)
  E <- pathway_order(gene_diff = D, scoreType = scoreType, eps = eps)
  return(E)
}

pathway_list <- c("Wnt signaling pathway", "Sonic Hedgehog (Shh) pathway", "Beta-catenin phosphorylation cascade", "Notch signaling pathway",
                  "FGF signaling pathway", "BMP receptor signaling", "BMP signaling and regulation", "BMP signaling pathway in stem cells",
                  "Retinoic acid receptor-mediated signaling", "Retinol metabolism")

## Evaluation method
evaluation_fun <- function(E) {
  res_mat <- matrix(NA, nrow = length(pathway_list), ncol = 3)
  rownames(res_mat) <- pathway_list
  colnames(res_mat) <- c("p-value", "p-adj", "position")
  for (iter in seq_len(length(pathway_list))) {
    index <- which(E$pathway == pathway_list[iter])
    res_mat[iter, 3] <- index
    res_mat[iter, 1:2] <- as.numeric(E[index, 2:3])
  }
  return(res_mat)
}

## with t-test filter
beta_use <- matrix(0, nrow = nrow(beta_mat), ncol = ncol(beta_mat)) 
beta_use[which(t_mat < 0.05)] <- beta_mat[which(t_mat < 0.05)] 
rownames(beta_use) <- rownames(beta_mat) 
colnames(beta_use) <- colnames(beta_mat)

## Manifold alignment
nGenes <- nrow(dta)
gene_name <- rownames(dta)
for (i in 2:30){
  mA_tmp <- res$res_large$mA[, 1:i]
  mA_X <- mA_tmp[1:nGenes, ]
  mA_Y <- mA_tmp[-(1:nGenes), ]
  ## return order of gene
  # gene_diff <- rowSums(abs((mA_X - mA_Y)))
  gene_diff <- rowSums((mA_X - mA_Y)^2)
  names(gene_diff) <- gene_name
  E <- pathway_order(gene_diff)
  print(length(which(E[, 3] < 0.05)))
  print(evaluation_fun(E))
}