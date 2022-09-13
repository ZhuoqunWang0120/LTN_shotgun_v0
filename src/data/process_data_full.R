#!/usr/bin/env Rscript

library(phyloseq)
library(ape)
library(data.tree)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
ntop = as.numeric(argv[2])
source(paste0(WORK_DIR,'/src/LTN/utils.R'))
system(paste0('mkdir -p ',WORK_DIR,'/cache'))
if (!file.exists(paste0(WORK_DIR,"/cache/ps_otu_tree.rds"))){
  raw_data = read.csv(paste0(WORK_DIR,"/res/diabimmune_t1d_16s_otu_table.txt"), sep = "\t", row.names = 1 )
  otutab0=apply(as.matrix(raw_data[,-ncol(raw_data)]),1:2,as.numeric)
  rownames(otutab0)=rownames(raw_data)
  colnames(otutab0)=colnames(raw_data)[-ncol(raw_data)]
  taxtab0=do.call(rbind,sapply(as.character(raw_data[,ncol(raw_data)]),function(x){strsplit(x,';  ')}))
  rownames(taxtab0)=rownames(otutab0)
  colnames(taxtab0)=c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
  tree=tax_tree(tax_table(taxtab0)) # slowest step, takes ~6hrs
  otu_preorder=tree$tip.label
  otutab0 = otu_table(t(otutab0[otu_preorder,]),taxa_are_rows = F)
  taxtab0 = tax_table(taxtab0[otu_preorder,])
  ps0=phyloseq(otutab0,taxtab0, tree)
  saveRDS(ps0,paste0(WORK_DIR,"/cache/ps_otu_tree.rds"))
}

# add sample_data
if (!file.exists(paste0(WORK_DIR,"/cache/ps_otu_full.rds"))){
  ps0 = readRDS(paste0(WORK_DIR,"/cache/ps_otu_tree.rds"))
  # sample data
  load(paste0(WORK_DIR,"/res/diabimmune_t1d_16s_metadata.rdata")) # md_16S
  rownames(md_16S)=as.character(md_16S$G_id)
  md_16S=md_16S[rownames(otu_table(ps0)),]
  ps0=merge_phyloseq(ps0,sample_data(md_16S))
  # subject data
  mmc2=readxl::read_excel(paste0(WORK_DIR,"/res/mmc2.xlsx"))
  mmc2=data.frame(mmc2)
  rownames(mmc2)=mmc2$ID..E.Espoo..Finland..T.Tartu..Estonia
  sampdat=sample_data(ps0)
  sampdat$sero_age=mmc2[as.character(sampdat$Subject_ID),"Age.at.sero.conversion..yrs"]*365
  sampdat$sero_age[is.na(sampdat$sero_age)]=max(sampdat$Age_at_Collection)+1
  sampdat$post_seroconversion=sampdat$Age_at_Collection>sampdat$sero_age
  sample_data(ps0)=sampdat
  saveRDS(ps0, paste0(WORK_DIR,"/cache/ps_otu_full.rds"))
}

ps0 = readRDS(paste0(WORK_DIR,"/cache/ps_otu_full.rds"))

# pruning
ra = t(apply(otu_table(ps0),1,function(x){x/sum(x)})) # nsample * nOTU
otu_top = names(sort(colSums(ra),decreasing = T)[1:ntop])
ps_top = prune_taxa(otu_top,ps0)
# reorder OTUs
otu_table(ps_top) = otu_table(ps_top)[,phy_tree(ps_top)$tip.label]
tax_table(ps_top) = tax_table(ps_top)[phy_tree(ps_top)$tip.label,]
saveRDS(ps_top,paste0(WORK_DIR,"/cache/ps_otu_",ntop,".rds"))
# yyl
yyl=seqtab2y(otu_table(ps_top),phy_tree(ps_top))
saveRDS(yyl,paste0(WORK_DIR,"/cache/yyl_otu_",ntop,".rds"))

