#!/usr/bin/env Rscript

rm(list=ls())
WORK_DIR = argv[1]
library(ggplot2)
library(dplyr)
library(reshape2)
RESULT_DIR = paste0(WORK_DIR,'/cache/tree/')
est_left = grep('est_para_left_lambda', list.files(RESULT_DIR, full.names = T), value = T)
for (f in est_left){
  system(paste0('mv ', f, ' ', RESULT_DIR, '/map2T1'))
}

K = 64
mu = rep(2, K - 1)
sigma = diag(K - 1)
est_mu = NULL
est_corr = NULL
files = NULL
for (lambda in c(0.1, 1, 3, 10, 30, 100)){
  rightfiles = paste0(RESULT_DIR, 'map2T1/est_right2left_lambda', lambda, 'K64_', 1:100, '.rds')
  balancedfiles = paste0(RESULT_DIR, 'map2T1/est_balanced2left_lambda', lambda, 'K64_', 1:100, '.rds')
  leftfiles = paste0(RESULT_DIR, 'map2T1/est_para_left_lambda', lambda, 'K64_', 1:100, '.rds')
  files = c(files, leftfiles, balancedfiles, rightfiles)
}

tictoc::tic()
for (i in seq_along(files)){
  if (file.exists(files[i])){
    est = readRDS(files[i])
    df1 = data.frame(mu = est[[1]])
    df2 = data.frame(corr = cov2cor(est[[2]])[upper.tri(est[[2]])])
    lambda = as.numeric(strsplit(strsplit(files[i],'lambda')[[1]][2],'K')[[1]][1])
    if (length(grep('right',strsplit(files[i],'_')[[1]]))>0){
      tree_structure = 'T3 (opposite)'
    } 
    else if (length(grep('balanced',strsplit(files[i],'_')[[1]]))>0){
      tree_structure = 'T2 (balanced)'
    }
    else{
      tree_structure = 'T1 (correct)'
    }
    df1$tree = tree_structure
    df2$tree = tree_structure
    df1$lambda = lambda
    df2$lambda = lambda
    est_mu = rbind(est_mu,df1)
    est_corr = rbind(est_corr,df2)
  }
}

est_mu$node_T1 = (1:63)
pair_index = matrix(0, K-1, K-1)
for (i in 1:(K-1)){
  for (j in 1:(K-1)){
    pair_index[i,j] = paste0(i-1,',',j-1)
  }
}
est_corr$pair_T1 = pair_index[upper.tri(pair_index)]
pairsplit = strsplit(est_corr$pair_T1,',')
est_corr$node1 = unlist(lapply(pairsplit, function(x){as.numeric(x[[1]])}))
est_corr$node2 = unlist(lapply(pairsplit, function(x){as.numeric(x[[2]])}))
tictoc::toc()

lam = 10

# boxplot
png(paste0(WORK_DIR, '/results/tree/boxplot_mu100.png'),width = 1500, height = 500)
est_mu %>% subset(lambda == lam) %>% 
  mutate(depth = as.factor(node_T1 - 1)) %>%
  mutate(tree = substr(tree,1,2)) %>%
  ggplot() +
  geom_boxplot(aes(x = depth, y = mu, color = tree)) +
  geom_abline(intercept = 2, slope = 0) +
  theme_bw()
dev.off()

png(paste0(WORK_DIR, '/results/tree/boxplot_corr100.png'),width = 1500, height = 500)
est_corr %>% subset(lambda == lam & ((node1 %in% 0:9 & node2 %in% 0:9) | (node1 %in% 53:62 & node2 %in% 53:62))) %>% 
  mutate(deep = factor(ifelse(node2 > 50, 'deep', 'shallow'), levels = c('shallow', 'deep'))) %>%
  mutate(tree = substr(tree,1,2)) %>%
  ggplot() +
  geom_boxplot(aes(x = pair_T1, y = corr, color = tree)) +
  facet_wrap(~deep, scales = 'free_x') +
  geom_abline(intercept = 0, slope = 0) +
  xlab('node pairs') +
  ylab('marginal correlation') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

# mse: change this into table 
mse_mu = NULL
size_mu = NULL
for (lam in unique(est_mu$lambda)){
  mse_mu = rbind(mse_mu, c(mean((est_mu %>% subset(lambda == lam & tree == 'T1 (correct)') %>% select(mu) %>% unlist() - 2)^2),
    mean((est_mu %>% subset(lambda == lam & tree == 'T2 (balanced)') %>% select(mu) %>% unlist() - 2)^2),
  mean((est_mu %>% subset(lambda == lam & tree == 'T3 (opposite)') %>% select(mu)  %>% unlist()  - 2)^2), lam))
  size_mu = rbind(size_mu, c(length((est_mu %>% subset(lambda == lam & tree == 'T1 (correct)') %>% select(mu) %>% unlist() - 2)^2),
                             length((est_mu %>% subset(lambda == lam & tree == 'T2 (balanced)') %>% select(mu) %>% unlist() - 2)^2),
                             length((est_mu %>% subset(lambda == lam & tree == 'T3 (opposite)') %>% select(mu)  %>% unlist()  - 2)^2), lam))
}
colnames(mse_mu) = c('T1','T2','T3','lambda')
colnames(size_mu) = c('T1','T2','T3','lambda')
mse_mu = mse_mu[order(mse_mu[,'lambda']),]
size_mu = size_mu[order(size_mu[,'lambda']),]


mse_corr = NULL
for (lam in unique(est_mu$lambda)){
  mse_corr = rbind(mse_corr, c(mean((est_corr %>% subset(lambda == lam & tree == 'T1 (correct)') %>% select(corr) %>% unlist() - 0)^2),
                               mean((est_corr %>% subset(lambda == lam & tree == 'T2 (balanced)') %>% select(corr) %>% unlist() - 0)^2),
                               mean((est_corr %>% subset(lambda == lam & tree == 'T3 (opposite)') %>% select(corr)  %>% unlist()  - 0)^2), lam))
}
colnames(mse_corr) = c('T1','T2','T3','lambda')
mse_corr
# mean((est_mu %>% subset(lambda == lam & tree == 'T1 (correct)') %>% select(mu) %>% unlist() %>% plogis() - plogis(2))^2)
# mean((est_mu %>% subset(lambda == lam & tree == 'T2 (balanced)') %>% select(mu) %>% unlist() %>% plogis() - plogis(2))^2)
# mean((est_mu %>% subset(lambda == lam & tree == 'T3 (opposite)') %>% select(mu)  %>% unlist() %>% plogis() - plogis(2))^2)

plt1 <- mse_mu %>% data.frame() %>% melt(id.vars = 'lambda') %>% mutate(tree = variable) %>% ggplot() + 
  geom_point(aes(x = as.factor(lambda), y = value, color = tree)) +
  xlab('lambda') + 
  ylab('MSE') + 
  theme_bw() +
  ggtitle('mean')
plt2 <- mse_corr %>% data.frame() %>% melt(id.vars = 'lambda') %>% mutate(tree = variable) %>% ggplot()+
  geom_point(aes(x = as.factor(lambda), y = value, color = tree)) +
  xlab('lambda') + 
  ylab('MSE') + 
  theme_bw() +
  ggtitle('marginal correlation')
png(paste0(WORK_DIR,'/results/tree/mse_tree.png'),width = 1000, height = 500)
gridExtra::grid.arrange(plt1,plt2,nrow = 1)
dev.off()
