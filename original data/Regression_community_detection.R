library(Matrix)
library(ggplot2)
# library(RSpectra)
# library(apcluster)
# library(FNN)
# library(NMF)

# ## load data
# m1 = Matrix(h5read("data/A1.mat", name = "data"))
# m2  = Matrix(h5read("data/A2.mat", name = "data"))
# m3  = Matrix(h5read("data/A3.mat", name = "data"))
# m4 = Matrix(h5read("data/A4.mat", name = "data"))
# m5 = Matrix(h5read("data/A5.mat", name = "data"))
# m6 = Matrix(h5read("data/A6.mat", name = "data"))
# m7 = Matrix(h5read("data/A7.mat", name = "data"))
# m8 = Matrix(h5read("data/A8.mat", name = "data"))
# m9 = Matrix(h5read("data/A9.mat", name = "data"))
# m10 = Matrix(h5read("data/A10.mat", name = "data"))
# n_gene <- nrow(m1)
# 
# ## Get response Y
# tmp1 <- as.numeric(m1)
# tmp2 <- as.numeric(m2)
# tmp3 <- as.numeric(m3)
# tmp4 <- as.numeric(m4)
# tmp5 <- as.numeric(m5)
# tmp6 <- as.numeric(m6)
# tmp7 <- as.numeric(m7)
# tmp8 <- as.numeric(m8)
# tmp9 <- as.numeric(m9)
# tmp10 <- as.numeric(m10)
# 
# time1 <- Sys.time()
# Y <- tmp1
# Y <- rbind(Y, tmp2)
# Y <- rbind(Y, tmp3)
# Y <- rbind(Y, tmp4)
# Y <- rbind(Y, tmp5)
# Y <- rbind(Y, tmp6)
# Y <- rbind(Y, tmp7)
# Y <- rbind(Y, tmp8)
# Y <- rbind(Y, tmp9)
# Y <- rbind(Y, tmp10)
# 
# 
# rm(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
# rm(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10)
# Y <- as(Y, "sparseMatrix")
# print(Sys.time() - time1) # Time difference of 2.475589 mins
# 
# saveRDS(Y, "data/largemat.rds")

m1 = readRDS("data/A1.rds")
n_gene <- nrow(m1)
gene_name <- read.table("data/genelist.txt")
gene_name <- gene_name$V1
Y <- readRDS("data/largemat.rds")
rm(m1)

source("code/utility.R")

##################################
#### Y = 1 + X + X^2   ###########
##################################

## Do regression
X <- matrix(c(rep(1, 10), (1:10), (1:10)^2), nrow = 10)
res <- my_regression(X, Y)
beta_coef <- res$coef
lm_t_p <- res$lm_t_p

# # use lm function to check results
# lm.fit <- lm(Y[, 7] ~ -1 + X)
# summary(lm.fit)
# beta_coef[, 7]
# lm_t_p[ , 7]
# lm_F_p[7]

## Focus on qudratic term first
beta_mat_quad <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[3, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p[2, ] < 0.05)
beta_mat_quad[index_p] <- beta_mat_tmp[index_p]

## Do community detection
index2commun <- function(index, gene_name, thred = 10){
  index_table <- which(table(index) > thred)
  index_commun <- list()
  for (i in 1:length(index_table)){
    index_commun[[i]] <- gene_name[which(index == as.numeric(names(index_table[i])))]
  }
  return(index_commun)
}

index_c <- which(colSums(abs(beta_mat_quad)) == 0)
index_r <- which(rowSums(abs(beta_mat_quad)) == 0)
index0 <- union(index_c, index_r)
index_use <- (1:n_gene)[-index0]
beta_mat_quad_use <- abs(beta_mat_quad[index_use, index_use])
gene_use <- gene_name[index_use]

## APC results.
# APC_res <- APC_fun(beta_mat_quad_use)
# print("Finish APC res.")
# K <- length(unique(APC_res))
# saveRDS(APC_res, "APC_quad.rds")
APC_res <- readRDS("results/cluster/APC_quad.rds")
APC_commun <- index2commun(index = APC_res, gene_name = gene_use, thred = 10)
saveRDS(APC_commun, "results/community/APC_commun_quad.rds")

## SCORE results.
# SCORE_res <- SCORE_fun(net = beta_mat_quad_use, K = K)
# print("Finish SCORE res.")
# saveRDS(SCORE_res, "SCORE_quad.rds")
SCORE_res <- readRDS("results/cluster/SCORE_quad.rds")
SCORE_commun <- index2commun(index = SCORE_res, gene_name = gene_use, thred = 10)
saveRDS(SCORE_commun, "results/community/SCORE_commun_quad.rds")

## LPA results.
# LPA_res <- LPA_fun(beta_mat_quad_use, K = 50)
# print("Finish LPA res.")
# saveRDS(LPA_res, "LPA_quad.rds")
LPA_res <- readRDS("results/cluster/LPA_quad.rds")
LPA_res <- LPA_res$community
LPA_commun <- index2commun(index = LPA_res, gene_name = gene_use, thred = 10)
saveRDS(LPA_commun, "results/community/LPA_commun_quad.rds")

# ## NMF results.
# NMF_res <- NMF_fun(beta_mat_quad_use, K = K)
# print("Finish NMF res.")
# saveRDS(NMF_res, "NMF_quad.rds")

## Calculate the center
beta2_tmp <- beta_coef[3, index_p]
beta1_tmp <- beta_coef[2, index_p]
center_quad <- -beta1_tmp / (2 * beta2_tmp)
summary(center_quad)
ggplot(data = NULL, aes(x = center_quad)) +
  geom_histogram(binwidth = 0.5)

# # function to split center
# index_to_community <- function(center_quad, low_center, up_center, index0, beta_mat_quad){
#   index_low <- which(center_quad < up_center)
#   index_high <- which(center_quad > low_center)
#   index_center <- intersect(index_low, index_high)
# 
#   index_use <- intersect((1:n_gene)[-index0], index_center)
#   beta_mat_quad_use <- abs(beta_mat_quad[index_use, index_use])
# 
#   res <- list()
#   res$index <- index_use
#   APC_res <- APC_fun(beta_mat_quad_use)
#   print("Finish APC res.")
# 
#   K <- length(unique(APC_res))
#   res$APC_res <- APC_res
# 
#   SCORE_res <- SCORE_fun(net = beta_mat_quad_use, K = K)
#   print("Finish SCORE res.")
#   res$SCORE_res <- SCORE_res
# 
#   LPA_res <- LPA_fun(beta_mat_quad_use, K = 25)
#   print("Finish LPA res.")
#   res$LPA_res <- LPA_res
#   return(res)
# }
# 
# res <- index_to_community(center_quad = center_quad, low_center = 2, up_center = 3, index0 = index0, beta_mat_quad = beta_mat_quad)
# saveRDS(res, "res2_3.rds")
# res <- index_to_community(center_quad = center_quad, low_center = 4, up_center = 5, index0 = index0, beta_mat_quad = beta_mat_quad)
# saveRDS(res, "res4_5.rds")
# res <- index_to_community(center_quad = center_quad, low_center = 6, up_center = 7, index0 = index0, beta_mat_quad = beta_mat_quad)
# saveRDS(res, "res6_7.rds")
# res <- index_to_community(center_quad = center_quad, low_center = 8, up_center = 9, index0 = index0, beta_mat_quad = beta_mat_quad)
# saveRDS(res, "res8_9.rds")

# get community
res2_3 <-readRDS("results/cluster/res2_3.rds")
res4_5 <-readRDS("results/cluster/res4_5.rds")
res6_7 <-readRDS("results/cluster/res6_7.rds")
res8_9 <-readRDS("results/cluster/res8_9.rds")
# APC results.
APC_commun_23 <- index2commun(index = res2_3$APC_res, gene_name = gene_name[res2_3$index], thred = 10)
APC_commun_45 <- index2commun(index = res4_5$APC_res, gene_name = gene_name[res4_5$index], thred = 10)
APC_commun_67 <- index2commun(index = res6_7$APC_res, gene_name = gene_name[res6_7$index], thred = 10)
APC_commun_89 <- index2commun(index = res8_9$APC_res, gene_name = gene_name[res8_9$index], thred = 10)
APC_res_split <- c(APC_commun_23, APC_commun_45, APC_commun_67, APC_commun_89)
saveRDS(APC_res_split, "results/community/APC_commun_quad_split.rds")
# SCORE results.
SCORE_commun_23 <- index2commun(index = res2_3$SCORE_res, gene_name = gene_name[res2_3$index], thred = 10)
SCORE_commun_45 <- index2commun(index = res4_5$SCORE_res, gene_name = gene_name[res4_5$index], thred = 10)
SCORE_commun_67 <- index2commun(index = res6_7$SCORE_res, gene_name = gene_name[res6_7$index], thred = 10)
SCORE_commun_89 <- index2commun(index = res8_9$SCORE_res, gene_name = gene_name[res8_9$index], thred = 10)
SCORE_res_split <- c(SCORE_commun_23, SCORE_commun_45, SCORE_commun_67, SCORE_commun_89)
saveRDS(SCORE_res_split, "results/community/SCORE_commun_quad_split.rds")
# LPA results.
LPA_commun_23 <- index2commun(index = res2_3$LPA_res$community, gene_name = gene_name[res2_3$index], thred = 10)
LPA_commun_45 <- index2commun(index = res4_5$LPA_res$community, gene_name = gene_name[res4_5$index], thred = 10)
LPA_commun_67 <- index2commun(index = res6_7$LPA_res$community, gene_name = gene_name[res6_7$index], thred = 10)
LPA_commun_89 <- index2commun(index = res8_9$LPA_res$community, gene_name = gene_name[res8_9$index], thred = 10)
LPA_res_split <- c(LPA_commun_23, LPA_commun_45, LPA_commun_67, LPA_commun_89)
saveRDS(LPA_res_split, "results/community/LPA_commun_quad_split.rds")

## Next is to consider linear term after removing quadratic term
index_lp <- setdiff(which(lm_t_p[1, ] < 0.05), index_p)
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[2, ], nrow = n_gene, ncol = n_gene)
beta_mat_linear[index_lp] <- beta_mat_tmp[index_lp]

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

# linear siginificant beta matrix
beta_mat_linear <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[2, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p < 0.05)
beta_mat_linear[index_p] <- beta_mat_tmp[index_p]

# Focus on the linear first
index_c <- which(colSums(abs(beta_mat_linear)) == 0)
index_r <- which(rowSums(abs(beta_mat_linear)) == 0)
index0 <- union(index_c, index_r)
index_use <- (1:n_gene)[-index0]
beta_mat_linear_use <- abs(beta_mat_linear[index_use, index_use])
gene_use <- gene_name[index_use]

## Do cummunity detection.
# APC results.
# APC_res <- APC_fun(beta_mat_linear_use)
# K <- length(unique(APC_res))
# print("Finish APC res.")
# saveRDS(APC_res, "APC_linear.rds")
APC_res <- readRDS("results/cluster/APC_linear.rds")
APC_commun <- index2commun(index = APC_res, gene_name = gene_use, thred = 10)
saveRDS(APC_commun, "results/community/APC_commun_linear.rds")

# SCORE results.
# SCORE_res <- SCORE_fun(net = beta_mat_linear_use, K = K)
# print("Finish SCORE res.")
# saveRDS(SCORE_res, "SCORE_linear.rds")
SCORE_res <- readRDS("results/cluster/SCORE_linear.rds")
SCORE_commun <- index2commun(index = SCORE_res, gene_name = gene_use, thred = 10)
saveRDS(SCORE_commun, "results/community/SCORE_commun_linear.rds")

# LPA results.
# LPA_res <- LPA_fun(beta_mat_linear_use, K = 25)
# print("Finish LPA res.")
# saveRDS(LPA_res, "LPA_linear.rds")
LPA_res <- readRDS("results/cluster/LPA_linear.rds")
LPA_res <- LPA_res$community
LPA_commun <- index2commun(index = LPA_res, gene_name = gene_use, thred = 10)
saveRDS(LPA_commun, "results/community/LPA_commun_linear.rds")


# NMF results.
# NMF_res <- NMF_fun(beta_mat_linear_use, K = K)
# print("Finish NMF res.")
# saveRDS(NMF_res, "NMF_linear.rds")


