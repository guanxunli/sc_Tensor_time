library(Matrix)
library(ggplot2)
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

get_KEGG_dis <- function(tsne_dta, db) {
  tsnePositions <- tsne_dta
  tsnePositions <- t(t(tsnePositions)/apply(tsnePositions, 2, function(X){max(abs(X))}))
  rownames(tsnePositions) <- toupper(rownames(tsnePositions))
  avgDistance <- lapply(db, function(gSet){
    gSet <- gSet[gSet %in% rownames(tsnePositions)]
    if(length(gSet) > 2){
      gSet <- tsnePositions[gSet,]
      gSet <- mahalanobis(gSet, colMeans(gSet), cov(tsnePositions))
      mean(gSet)  
    }
  })
  return(avgDistance)
}
set.seed(1234)

## check benchmark distance
# m1 <- readRDS("data/A1.rds")
# m2 <- readRDS("data/A2.rds")
# m3 <- readRDS("data/A3.rds")
m4 <- readRDS("data/A4.rds")
# m5 <- readRDS("data/A5.rds")
# m6 <- readRDS("data/A6.rds")
m7 <- readRDS("data/A7.rds")
# m8 <- readRDS("data/A8.rds")
# m9 <- readRDS("data/A9.rds")
# m10 <- readRDS("data/A10.rds")
# mavg <- (m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10) / 10
# saveRDS(mavg, "data/Aavg.rds")
mavg <- readRDS("data/Aavg.rds")

m4_tsne <- my_tsne(m4, gene_list = gene_name)
m7_tsne <- my_tsne(m7, gene_list = gene_name)
mavg_tsne <- my_tsne(mavg, gene_list = gene_name)
rm(m4, m7, mavg)
m4_KEGG <- get_KEGG_dis(m4_tsne, db = KEGG)
m7_KEGG <- get_KEGG_dis(m7_tsne, db = KEGG)
mavg_KEGG <- get_KEGG_dis(mavg_tsne, db = KEGG)

KEGG_figure <- function(KEGG_mat) {
    
  gSets <- intersect(names(m4_KEGG), names(m7_KEGG))
  gSets <- intersect(gSets, names(mavg_KEGG))
  gSets <- intersect(gSets, names(KEGG_mat))
  
  m4Dist <- unlist(m4_KEGG[gSets])
  m7Dist <- unlist(m7_KEGG[gSets])
  mavgDist <- unlist(mavg_KEGG[gSets])
  betaDist <- unlist(KEGG_mat[gSets])
  
  m4Dist <- data.frame(dist = m4Dist, method = 'm4')
  m7Dist <- data.frame(dist = m7Dist, method = 'm7')
  mavgDist <- data.frame(dist = mavgDist, method = 'mavg')
  betaDist <- data.frame(dist = betaDist, method = 'beta')

  allDist <- rbind(m4Dist, m7Dist, mavgDist, betaDist)
  allDist$method <- factor(allDist$method, levels = c('m4', 'm7', 'mavg', 'beta'))

  kTest <- kruskal.test(allDist$dist, allDist$method)
  kTest <- as.numeric(formatC(kTest$p.value, digits = 3))

  A <- ggplot(allDist, aes(dist, color = method)) + stat_ecdf() + 
    theme_bw() + 
    xlab(Avg~(Mahalinobis~Distance)) + 
    ylab('Empirical Cumulative\n Density Function (ECDF)') + 
    labs(title = paste0('Kruskal test P-value = ', kTest))

  B <- ggplot(allDist, aes(method, dist, color = method)) + geom_boxplot(outlier.color = NA) + 
    geom_jitter(alpha = 0.1, cex = 0.2) + 
    theme_bw() + 
    xlab('Method') + 
    ylab(Avg~(Mahalinobis~Distance)) +
    labs(title = paste0('Kruskal test P-value = ', kTest))

  res <- list()
  res$A <- A
  res$B <- B
  return(res)
}

res_KEGG_list <- list()
res_KEGG_factor <- c("m4", "m7", "mavg")
res_KEGG_list[[1]] <- m4_KEGG
res_KEGG_list[[2]] <- m7_KEGG
res_KEGG_list[[3]] <- mavg_KEGG
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
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

# res_beta_quad <- KEGG_figure(beta_KEGG)
# png('beta_quad_withoutfilter.png', width = 1600, height = 1000, res = 300)
# res_beta_quad$B
# dev.off()

res_KEGG_list[[4]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "quad_withoutp")

# Take abs
beta_mat_quad <- abs(beta_mat_quad)

beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

res_KEGG_list[[5]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "quad_withoutp_abs")

## With p-value filtering
beta_mat_quad <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[3, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p[2, ] < 0.05)
beta_mat_quad[index_p] <- beta_mat_tmp[index_p]

beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

# res_beta_quad <- KEGG_figure(beta_KEGG)
# png('beta_quad_withfilter.png', width = 1600, height = 1000, res = 300)
# res_beta_quad$B
# dev.off()

res_KEGG_list[[6]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "quad_with_p")

# Take abs
beta_mat_quad <- abs(beta_mat_quad)
beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)
res_KEGG_list[[7]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "quad_with_p_abs")

## Next is to combine linear term after removing quadratic term
index_lp <- setdiff(which(lm_t_p[1, ] < 0.05), index_p)
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[2, ], nrow = n_gene, ncol = n_gene)
beta_mat_linear[index_lp] <- beta_mat_tmp[index_lp]

index_c <- which(colSums(abs(beta_mat_linear)) == 0)
index_r <- which(rowSums(abs(beta_mat_linear)) == 0)
index0 <- union(index_c, index_r)

index_cq <- which(colSums(abs(beta_mat_quad)) == 0)
index_rq <- which(rowSums(abs(beta_mat_quad)) == 0)
index0q <- union(index_cq, index_rq)

index0 <- union(index0, index0q)
index_use <- (1:n_gene)[-index0]
gene_use <- gene_list[index_use]

beta_mat <- cbind(beta_mat_quad[index_use, index_use], beta_mat_linear[index_use, index_use])
beta_mat <- as.matrix(beta_mat)

# Do t-sne
mat_svd <- svds(beta_mat, k = 100)
mat_pca <- beta_mat %*% mat_svd$v
mat_tsne <- Rtsne(mat_pca, dim = 2, initial_dims = 100, pca = FALSE, verbose = FALSE)
beta_tsne <- mat_tsne$Y
rownames(beta_tsne) <- gene_use

beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

res_beta_quad <- KEGG_figure(beta_KEGG)
png('beta_quad_linear_withfilter.png', width = 1600, height = 1000, res = 300)
res_beta_quad$B
dev.off()

## consider 0 - 1 case
beta_mat_quad <- matrix(0, nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p[2, ] < 0.05)
beta_mat_quad[index_p] <- 1

beta_tsne <- my_tsne(mat = beta_mat_quad, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

# res_beta_quad <- KEGG_figure(beta_KEGG)
# png('beta_quad_01.png', width = 1600, height = 1000, res = 300)
# res_beta_quad$B
# dev.off()

res_KEGG_list[[8]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "quad_01")

index_lp <- setdiff(which(lm_t_p[1, ] < 0.05), index_p)
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_linear[index_lp] <- 1

beta_mat <- beta_mat_quad + beta_mat_linear
beta_tsne <- my_tsne(mat = beta_mat, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

# res_beta_quad <- KEGG_figure(beta_KEGG)
# png('beta_quad_linear_01.png', width = 1600, height = 1000, res = 300)
# res_beta_quad$B
# dev.off()

res_KEGG_list[[9]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "quad_linear_01")


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
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

# res_beta_quad <- KEGG_figure(beta_KEGG)
# png('beta_linear_withoutfilter.png', width = 1600, height = 1000, res = 300)
# res_beta_quad$B
# dev.off()

res_KEGG_list[[10]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "linear_without_p")

# Take abs
beta_mat_linear <- abs(beta_mat_linear)

beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

res_KEGG_list[[11]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "linear_without_p_abs")

## linear siginificant beta matrix
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[2, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p < 0.05)
beta_mat_linear[index_p] <- beta_mat_tmp[index_p]

beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

# res_beta_quad <- KEGG_figure(beta_KEGG)
# png('beta_linear_withfilter.png', width = 1600, height = 1000, res = 300)
# res_beta_quad$B
# dev.off()

res_KEGG_list[[12]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "linear_with_p")

# Take abs
beta_mat_linear <- abs(beta_mat_linear)

beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

res_KEGG_list[[13]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "linear_with_p_abs")

# 0 - 1 case
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p < 0.05)
beta_mat_linear[index_p] <- 1

beta_tsne <- my_tsne(mat = beta_mat_linear, gene_list = gene_name)
beta_KEGG <- get_KEGG_dis(tsne_dta = beta_tsne, db = KEGG)

# res_beta_quad <- KEGG_figure(beta_KEGG)
# png('beta_linear_01.png', width = 1600, height = 1000, res = 300)
# res_beta_quad$B
# dev.off()

res_KEGG_list[[14]] <- beta_KEGG
res_KEGG_factor <- c(res_KEGG_factor, "linear_01")

#### Combine all results!
gSets <- names(res_KEGG_list[[1]])
for (i in 2:length(res_KEGG_list)){
  gSets <- intersect(gSets, names(res_KEGG_list[[i]]))
}

m4Dist <- unlist(res_KEGG_list[[1]][gSets])
m4Dist <- data.frame(dist = m4Dist, method = 'm4')
allDist <- m4Dist
for (i in 2:length(res_KEGG_list)){
  tmp_Dist <- unlist(res_KEGG_list[[i]][gSets])
  tmp_Dist <- data.frame(dist = tmp_Dist, method = res_KEGG_factor[i])
  allDist <- rbind(allDist, tmp_Dist)
}

allDist$method <- factor(allDist$method, levels = res_KEGG_factor)

kTest <- kruskal.test(allDist$dist, allDist$method)
kTest <- as.numeric(formatC(kTest$p.value, digits = 3))

ggplot(allDist, aes(dist, color = method)) + stat_ecdf() + 
    theme_bw() + 
    xlab(Avg~(Mahalinobis~Distance)) + 
    ylab('Empirical Cumulative\n Density Function (ECDF)') + 
    labs(title = paste0('Kruskal test P-value = ', kTest))

png('all_beta.png', width = 2500, height = 1500, res = 300)
ggplot(allDist, aes(method, dist, color = method)) + geom_boxplot(outlier.color = NA) + 
    geom_jitter(alpha = 0.1, cex = 0.2) + 
    theme_bw() + 
    xlab('Method') + 
    ylab(Avg~(Mahalinobis~Distance)) +
    labs(title = paste0('Kruskal test P-value = ', kTest)) +
    theme(axis.text.x = element_text(angle = -90, hjust=0))
dev.off()
