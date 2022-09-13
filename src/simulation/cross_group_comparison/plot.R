#!/usr/bin/env Rscript
library(phyloseq)
library(ape)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
source(paste0(WORK_DIR,'/src/LTN/utils.R'))
K=100
png(paste0(WORK_DIR,'/results/cross_group_comparison/ROC.png'),height = 500,width = 1000)
par(mfrow=c(1,2),pty='s')
# colPal= c(2,1,3,4,6)
colPal = 1:3
for (i in 1:2){
  s=c('single_otu','multi_otu')[i]
  resultdir=paste0(WORK_DIR,'/results/cross_group_comparison/',s)
  # roclist=sapply(c('r2.rds','r5.rds','sparse_lambda10.rds','diagonal_'),function(x){readRDS(grep(x,list.files(resultdir,full.names = T),value = T))})
  roclist=sapply(c('sparse_lambda10.rds','r2.rds','r5.rds'),function(x){
    filenam = grep(x,list.files(resultdir,full.names = T),value = T)
    if (length(filenam) == 1){
      return(readRDS(filenam))
    }
  })
  #for (j in 1:5){
  for (j in 1:3){
    if (!is.null(roclist[[j]])){
      plot(roclist[[j]],add=(j!=1),col=colPal[j],main=paste0('scenario ',i))
    }
  }
  legend('bottomright',
         # legend = c('LTN (lambda=10)','LTN (lambda=0.1)','LTN (diagonal)','DirFactor (r=2)','DirFactor (r=5)'),
         legend = c('LTN (lambda=10)', 'DirFactor (r=2)','DirFactor (r=5)'),
         fill=colPal)
}
dev.off()

