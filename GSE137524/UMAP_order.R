pathway_order <- function(dta, scoreType = "pos", eps = 0){
  BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
  set.seed(1)
  E <- fgseaMultilevel(BIOP, dta[!grepl('^RPL|^RPS|^RP[[:digit:]]|^MT-',names(dta))],
                       scoreType = scoreType, eps = eps)
  E <- E[order(E$pval, decreasing = FALSE),]
  E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
  index <- which(E$pathway == "EGFR1 pathway")
  res <- matrix(NA, nrow = 1, ncol = 3)
  tmp <- E[index, 1:3]
  rownames(res) <- tmp[1, 1]
  colnames(res) <- c("p-value", "p-adj", "position")
  res[1,3] <- index
  res[1,1:2] <- as.numeric(tmp[, 2:3])
  return(res)
}

UMAP_order <- function(dta, scoreType = "pos", eps = 0){
  set.seed(1)
  dta_svd <- RSpectra::svds(dta, k = 50)
  dta_pcscore <-  t(t(dta_svd$u) * dta_svd$d)
  UMAP <- umap(dta_pcscore)
  rownames(UMAP) <- rownames(dta)
  ## GSEA
  mDistance <- mahalanobis(UMAP, colMeans(UMAP), cov(UMAP))
  res <- pathway_order(mDistance, scoreType = scoreType, eps = eps)
  return(res)
}