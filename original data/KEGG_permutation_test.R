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

# my_tsne <- function(mat, gene_list) {
#   # remove 0 row and columns
#   mat <- as.matrix(mat)
#   index_c <- which(colSums(abs(mat)) == 0)
#   index_r <- which(rowSums(abs(mat)) == 0)
#   index0 <- union(index_c, index_r)
#   n_gene <- length(gene_list)
#   index_use <- (1:n_gene)[-index0]
#   gene_use <- gene_list[index_use]
#   mat_use <- mat[index_use, index_use]
#   mat_use <- exp(-abs(mat_use))
#   diag(mat_use) <- 0
#   
#   # Do t-sne
#   index <- matrix(1, nrow = nrow(mat_use), ncol = ncol(mat_use))
#   mat_tsne <- Rtsne_neighbors(index, distance = mat_use, dim = 2, verbose = FALSE, check_duplicates = FALSE)
#   mat_tsne <- mat_tsne$Y
#   rownames(mat_tsne) <- gene_use
#   return(mat_tsne)
# }

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
  KEGG_res <- KEGG_res[order(KEGG_res$p_value, decreasing = FALSE), ]
  return(KEGG_res)
}

set.seed(1234)
## check benchmark distance
mavg <- readRDS("data/Aavg.rds")
mavg_tsne <- my_tsne(mavg, gene_list = gene_name)
rm(mavg)
mavg_KEGG <- KEGG_permutaion_test(mavg_tsne, db = KEGG)
write.table(mavg_KEGG, "mavg_KEGG.txt")

##################################
#### Y = 1 + X + X^2   ###########
##################################

## Do regression
X <- matrix(c(rep(1, 10), (1:10), (1:10)^2), nrow = 10)
res <- my_regression(X, Y)
beta_coef <- res$coef
lm_t_p <- res$lm_t_p

#### Focus on qudratic term first
## without filter
beta_mat_quad <- matrix(beta_coef[3, ], nrow = n_gene, ncol = n_gene)

beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "quad_without_p.txt")

# Take abs
beta_mat_quad <- abs(beta_mat_quad)

beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "quad_without_p_abs.txt")

## With p-value filtering
beta_mat_quad <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[3, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p[2, ] < 0.05)
beta_mat_quad[index_p] <- beta_mat_tmp[index_p]

beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "quad_with_p.txt")

# Take abs
beta_mat_quad <- abs(beta_mat_quad)
beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "quad_with_p_abs.txt")

#####################################
########## Y = 1 + X  ###############
#####################################

# Do regress
X <- matrix(c(rep(1, 10), (1:10)), nrow = 10)
res <- my_regression(X = X, Y = Y)
beta_coef <- res$coef
lm_t_p <- res$lm_t_p

# # lm function to check results
# lm.fit <- lm(Y[, 7] ~ -1 + X)
# summary(lm.fit)
# beta_coef[, 7]
# lm_t_p[, 7]

## Without filter
beta_mat_linear <- matrix(beta_coef[2, ], nrow = n_gene, ncol = n_gene)

beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "linear_without_p.txt")

# Take abs
beta_mat_linear <- abs(beta_mat_linear)
beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "linear_without_p_abs.txt")

## linear siginificant beta matrix
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[2, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p < 0.05)
beta_mat_linear[index_p] <- beta_mat_tmp[index_p]

beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "linear_with_p.txt")

# Take abs
beta_mat_linear <- abs(beta_mat_linear)
beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- KEGG_permutaion_test(tsne_dta = beta_tsne, db = KEGG)
write.table(beta_KEGG, "linear_with_p_abs.txt")

############################### Check results #####################################
mavg_KEGG <- read.table("results/KEGG_test/pathway/original/mavg_KEGG.txt", header = TRUE)
quad_without_p <- read.table("results/KEGG_test/pathway/original/quad_without_p.txt", header = TRUE)
quad_without_p_abs <- read.table("results/KEGG_test/pathway/original/quad_without_p_abs.txt", header = TRUE)
quad_with_p <- read.table("results/KEGG_test/pathway/original/quad_with_p.txt", header = TRUE)
quad_with_p_abs <- read.table("results/KEGG_test/pathway/original/quad_with_p_abs.txt", header = TRUE)
linear_without_p <- read.table("results/KEGG_test/pathway/original/linear_without_p.txt", header = TRUE)
linear_without_p_abs <- read.table("results/KEGG_test/pathway/original/linear_without_p_abs.txt", header = TRUE)
linear_with_p <- read.table("results/KEGG_test/pathway/original/linear_with_p.txt", header = TRUE)
linear_with_p_abs <- read.table("results/KEGG_test/pathway/original/linear_with_p_abs.txt", header = TRUE)

KEGG_list <- list()
KEGG_list[[1]] <- mavg_KEGG
KEGG_list[[2]] <- quad_without_p
KEGG_list[[3]] <- quad_without_p_abs
KEGG_list[[4]] <- quad_with_p
KEGG_list[[5]] <- quad_with_p_abs
KEGG_list[[6]] <- linear_without_p
KEGG_list[[7]] <- linear_without_p_abs
KEGG_list[[8]] <- linear_with_p
KEGG_list[[9]] <- linear_with_p_abs

for (i in 1:9){
  tmp <- KEGG_list[[i]]
  tmp$p.adj <- p.adjust(tmp$p.value)
  KEGG_list[[i]] <- tmp
}

num_intersect <- function(dta1, dta2){
  
  tmp1 <- rownames(dta1)[which(dta1$p.adj < 0.05)]
  tmp2 <- rownames(dta2)[which(dta2$p.adj < 0.05)]
  
  # tmp1 <- rownames(dta1)[which(dta1$p.value < 0.05)]
  # tmp2 <- rownames(dta2)[which(dta2$p.value < 0.05)]
  
  # tmp1 <- rownames(dta1)[1:100]
  # tmp2 <- rownames(dta2)[1:100]
  
  return(length(intersect(tmp1, tmp2)))
}

num_mat <- matrix(NA, nrow = 9, ncol = 9)
for (i in 1:9){
  for (j in i:9){
    num_mat[i, j] <- num_intersect(KEGG_list[[i]], KEGG_list[[j]])
  }
}
num_mat

## more tests
pathway_mavg <- rownames(mavg_KEGG)[which(KEGG_list[[1]]$p.adj < 0.05)]
pathway_quad_with_p_abs <- rownames(quad_with_p_abs)[which(KEGG_list[[5]]$p.adj < 0.05)]
pathway_linear_with_p_abs <- rownames(linear_with_p_abs)[which(KEGG_list[[9]]$p.adj < 0.05)]

inter_quad_linear <- intersect(pathway_quad_with_p_abs, pathway_linear_with_p_abs)
setdiff(pathway_quad_with_p_abs, inter_quad_linear) %in% pathway_mavg
setdiff(pathway_linear_with_p_abs, inter_quad_linear) %in% pathway_mavg



