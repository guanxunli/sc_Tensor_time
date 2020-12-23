library(Matrix)
library(ggplot2)
library(patchwork)
library(gplots)

m1 = readRDS("data/A1.rds")
n_gene <- nrow(m1)
gene_name <- read.table("data/genelist.txt")
gene_name <- gene_name$V1
Y <- readRDS("data/largemat.rds")
rm(m1)

source("code/utility.R")

## Do regression
X <- matrix(c(rep(1, 10), (1:10), (1:10)^2), nrow = 10)
res <- my_regression(X, Y)
beta_coef <- res$coef
lm_t_p <- res$lm_t_p

## Focus on qudratic term first
beta_mat_quad <- matrix(0, nrow = n_gene, ncol = n_gene)
beta_mat_tmp <- matrix(beta_coef[3, ], nrow = n_gene, ncol = n_gene)
index_p <- which(lm_t_p[2, ] < 0.05)
beta_mat_quad[index_p] <- beta_mat_tmp[index_p]
rownames(beta_mat_quad) <- gene_name
colnames(beta_mat_quad) <- gene_name



LPA_commun <- readRDS("results/community/LPA_commun_quad_split.rds")
list_num <- 9
gene_name <- LPA_commun[[list_num]]
beta_mat <- beta_mat_quad[gene_name, gene_name]
png(paste0('results/criterion/heatmap_LPA_split', list_num, '.png'), width = 2500, height = 1000, res = 300, pointsize = 10)
heatmap.2(as.matrix(beta_mat), margins=c(5,8), key=TRUE, symkey=FALSE, density.info="none", trace="none")
dev.off()





