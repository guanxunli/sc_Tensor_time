BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')

pathway_order <- function(gene_diff, scoreType = "pos", eps = 0) {
  set.seed(1)
  names(gene_diff) <- toupper(names(gene_diff))
  E <- fgseaMultilevel(BIOP, gene_diff[!grepl('^RPL|^RPS|^RP[[:digit:]]|^MT-',names(gene_diff))],
                       scoreType = scoreType, eps = eps)
  E <- E[order(E$pval, decreasing = FALSE),]
  E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
  return(E)
}

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
  E <- pathway_order(gene_diff = mDistance, scoreType = scoreType, eps = eps)
  return(E)
}

pathway_list <- c("Wnt signaling pathway", "Sonic Hedgehog (Shh) pathway", "Beta-catenin phosphorylation cascade", "Notch signaling pathway",
                  "FGF signaling pathway", "BMP receptor signaling", "BMP signaling and regulation", "BMP signaling pathway in stem cells",
                  "Retinoic acid receptor-mediated signaling", "Retinol metabolism")

evaluation_fun <- function(E) {
  res_mat <- matrix(NA, nrow = length(pathway_list), ncol = 3)
  rownames(res_mat) <- pathway_list
  colnames(res_mat) <- c("p-value", "p-adj", "position")
  for (iter in seq_len(length(pathway_list))) {
    index <- which(E$pathway == pathway_list[iter])
    res_mat[iter, 3] <- index
    res_mat[iter, 1:2] <- as.numeric(E[index, 1:2])
  }
  return(res_mat)
}