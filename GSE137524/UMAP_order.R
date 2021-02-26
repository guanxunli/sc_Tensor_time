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