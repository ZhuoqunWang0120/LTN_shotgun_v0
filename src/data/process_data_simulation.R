#!/usr/bin/env Rscript

library(phyloseq)
library(data.tree)
library(ape)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
ntop = as.numeric(argv[2])
ps_top=readRDS(paste0(WORK_DIR,"/cache/ps_otu_",ntop,".rds"))
md_16S=sample_data(ps_top)
samples_finland_control=rownames(md_16S)[md_16S$Case_Control=='control' & md_16S$Country=='Finland']
ps_finland_control=prune_samples(samples_finland_control,ps_top)
tree=phy_tree(ps_finland_control)
cnt=matrix(otu_table(ps_finland_control),ncol=ncol(otu_table(ps_finland_control)),dimnames = dimnames(otu_table(ps_finland_control)))
sampdat=sample_data(ps_finland_control)
age=sampdat$Age_at_Collection
age_s=(age-mean(age))/sd(age)
Xadjust=matrix(age_s,ncol=1)
Xadjust=cbind(1,Xadjust)
grouplabel=as.numeric(sampdat$Subject_ID)
g=max(grouplabel)
# gprior_m=1
tree=as.Node(phy_tree(ps_finland_control))
tree$Set(name=(ntop+1):(2*ntop-1),traversal = 'pre-order',filterFun = isNotLeaf)
tree=as.phylo(tree)
# identical(tree$tip.label,colnames(otu_table(ps_top)))
ra_ntop=apply(cnt, 1, function(x){x/sum(x)})
top20_ra=names(sort(rowSums(ra_ntop),decreasing = T)[1:20])
saveRDS(list(cnt=cnt,Xadjust=Xadjust,grouplabel=grouplabel,g=g,tree=tree,top20_ra=top20_ra),paste0(WORK_DIR,"/cache/ps_sim_",ntop,".rds"))


