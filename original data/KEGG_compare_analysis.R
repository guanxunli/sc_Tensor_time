library(Matrix)
library(RSpectra)
library(Rtsne)

## load data
m1 = readRDS("data/A1.rds")
n_gene <- nrow(m1)
rm(m1)

gene_name <- read.table("data/genelist.txt")
gene_name <- gene_name$V1
Y <- readRDS("data/largemat.rds")
source("code/utility.R")

##################################
#### load needed #################
##################################
## KEGG test
library(fgsea)
KEGG <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')

my_tsne <- function(mat, gene_list) {
  # remove 0 row and columns
  mat <- as.matrix(mat)
  index_c <- which(colSums(abs(mat)) == 0)
  index_r <- which(rowSums(abs(mat)) == 0)
  index0 <- union(index_c, index_r)
  n_gene <- length(gene_list)
  index_use <- (1:n_gene)[-index0]
  gene_use <- gene_list[index_use]
  mat_use <- mat[index_use, index_use]
  
  # Do t-sne
  mat_svd <- svds(mat_use, k = 100)
  mat_pca <- mat_use %*% mat_svd$v
  mat_tsne <- Rtsne(mat_pca, dim = 2, initial_dims = 100, pca = FALSE, verbose = FALSE, check_duplicates = FALSE)
  mat_tsne <- mat_tsne$Y
  rownames(mat_tsne) <- gene_use
  return(mat_tsne)
}

KEGG_permutaion_test <- function(tsne_dta, db = KEGG, iter_per = 10000){
  tsnePositions <- tsne_dta
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  npathway <- length(db)
  KEGG_res <- matrix(NA, nrow = npathway, ncol = 3)
  rownames(KEGG_res) <- names(db)
  colnames(KEGG_res) <- c("average distance", "p-value", "number of gene")
  ngene <- nrow(tsne_dta)
  
  for (iter in 1:npathway){
    gSet <- db[[iter]]
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    nGene <- length(gSet)
    if (nGene > 2){
      KEGG_res[iter, 3] <- nGene
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      dDis <- mean(gSet)
      KEGG_res[iter, 1] <- dDis
      dDis_per <- rep(NA, iter_per)
      for (i in 1:iter_per){
        gSet <- tsnePositions[sample(1:ngene, nGene),]
        gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
        dDis_per[i] <- mean(gSet)
      }
      p_value <- length(which(dDis_per < dDis)) / iter_per
      KEGG_res[iter, 2] <- p_value
    }
  }
  KEGG_res <- na.omit(KEGG_res)
  KEGG_res <- KEGG_res[order(KEGG_res[, 2], decreasing = FALSE), ]
  return(KEGG_res)
}

set.seed(1234)
## check benchmark distance
mavg <- readRDS("data/Aavg.rds")
mavg_tsne <- my_tsne(mavg, gene_list = gene_name)
rm(mavg)
mavg_KEGG <- KEGG_permutaion_test(mavg_tsne, db = KEGG)
write.table(mavg_KEGG, "results/mavg_KEGG.txt")

##################################
#### Y = 1 + X + X^2   ###########
##################################

## Do regression
X <- matrix(c(rep(1, 10), (1:10), (1:10)^2), nrow = 10)
res <- my_regression(X, Y)
beta_coef <- res$coef
lm_t_p <- res$lm_t_p

#### Focus on qudratic term first
## With p-value filtering
beta_mat_quad <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[3, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p[2, ] < 0.05)
## With center filtering
beta2_tmp <- beta_coef[3, index_p]
beta1_tmp <- beta_coef[2, index_p]
center_quad <- -beta1_tmp / (2 * beta2_tmp)
index_l <- which(center_quad > 3)
index_r <- which(center_quad < 8)
index_p <- index_p[intersect(index_l, index_r)]
beta_mat_quad[index_p] <- beta_mat_tmp[index_p]

beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "quad_p.txt")

# Take abs
beta_mat_quad <- abs(beta_mat_quad)
beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "quad_p_abs.txt")

#### Next is to focus on the linear term
## With p-value filtering
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[2, ], nrow = n_gene, ncol = n_gene)
index_l <- setdiff(which(lm_t_p[1, ] < 0.05), index_p)
beta_mat_linear[index_l] <- beta_mat_tmp[index_l]

beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "linear_p.txt")

# Take abs
beta_mat_linear <- abs(beta_mat_linear)
beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "linear_p_abs.txt")

################################ Compare results ##############################
mavg_KEGG <- read.table("results/KEGG_test/pathway/compare/mavg_KEGG.txt", header = TRUE)
mavg_KEGG$p.adj <- p.adjust(mavg_KEGG$p.value)
mavg_pathway <- rownames(mavg_KEGG)[which(mavg_KEGG$p.adj < 0.05)]

## with out abs
quad_KEGG <- read.table("results/KEGG_test/pathway/compare/quad_p.txt", header = TRUE)
quad_KEGG$p.adj <- p.adjust(quad_KEGG$p.value)
quad_pathway <- rownames(quad_KEGG)[which(quad_KEGG$p.adj < 0.05)]

qlinear_KEGG <- read.table("results/KEGG_test/pathway/compare/linear_p.txt", header = TRUE)
qlinear_KEGG$p.adj <- p.adjust(qlinear_KEGG$p.value)
qlinear_pathway <- rownames(qlinear_KEGG)[which(qlinear_KEGG$p.adj < 0.05)]

print(c(length(quad_pathway), length(qlinear_pathway), length(intersect(quad_pathway, qlinear_pathway))))
quad_linear_pathway <- union(quad_pathway, qlinear_pathway)
print(c(length(quad_linear_pathway), length(mavg_pathway), length(intersect(quad_linear_pathway, mavg_pathway))))

linear_p <- read.table("results/KEGG_test/pathway/compare/linear_with_p.txt", header = TRUE)
linear_p$p.adj <- p.adjust(linear_p$p.value)
linear_pathway <- rownames(linear_p)[which(linear_p$p.adj < 0.05)]
print(c(length(quad_linear_pathway), length(linear_pathway), length(intersect(quad_linear_pathway, linear_pathway))))

## with abs
quad_KEGG <- read.table("results/KEGG_test/pathway/compare/quad_p_abs.txt", header = TRUE)
quad_KEGG$p.adj <- p.adjust(quad_KEGG$p.value)
quad_pathway <- rownames(quad_KEGG)[which(quad_KEGG$p.adj < 0.05)]

qlinear_KEGG <- read.table("results/KEGG_test/pathway/compare/linear_p_abs.txt", header = TRUE)
qlinear_KEGG$p.adj <- p.adjust(qlinear_KEGG$p.value)
qlinear_pathway <- rownames(qlinear_KEGG)[which(qlinear_KEGG$p.adj < 0.05)]

print(c(length(quad_pathway), length(qlinear_pathway), length(intersect(quad_pathway, qlinear_pathway))))
quad_linear_pathway <- union(quad_pathway, qlinear_pathway)
print(c(length(quad_linear_pathway), length(mavg_pathway), length(intersect(quad_linear_pathway, mavg_pathway))))

linear_p <- read.table("results/KEGG_test/pathway/compare/linear_with_p_abs.txt", header = TRUE)
linear_p$p.adj <- p.adjust(linear_p$p.value)
linear_pathway <- rownames(linear_p)[which(linear_p$p.adj < 0.05)]
print(c(length(quad_linear_pathway), length(linear_pathway), length(intersect(quad_linear_pathway, linear_pathway))))


