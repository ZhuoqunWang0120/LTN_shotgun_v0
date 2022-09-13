#!/usr/bin/env Rscript

library(reshape2)
library(ggplot2)
library(dplyr)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
nsim1=argv[2]
result_dir = paste0(WORK_DIR,"/results/covariance_estimation/")
files = list.files(result_dir,full.names = T)
LOSS = data.frame()
for (lambda in c(0.1,1,10)){
  for (mu_var in c(1,5,50,100)){
    loss = melt(readRDS(grep(paste0('loss_LTN_lambda_',lambda,'_mu_var',mu_var,'.rds'),files,value = T)))
    loss$lambda = lambda
    loss$mu_var = mu_var
    loss$model = 'LTN'
    LOSS = rbind(LOSS, loss)
  }
}
loss = melt(readRDS(grep('COAT',files,value = T)))
loss$model = 'COAT'
loss$lambda = NA
loss$mu_var = NA
loss$Var1 = sapply(loss$Var1,function(x){if(x == 'L_inf'){return('L_infinity')} else{return(as.character(x))}})
#
LOSS_coat = LOSS
LOSS_coat$value = loss$value
LOSS_coat$model = 'COAT'
LOSS_coat$lambda = -999
LOSS = rbind(LOSS,LOSS_coat)
LOSS = LOSS %>% group_by(Var1) %>% mutate(y_min = min(value), y_max = max(value) * 1.2) %>% ungroup()
LOSS$interaction = interaction(LOSS$model,LOSS$lambda)
LOSS$linetype = ifelse(LOSS$model == 'COAT', 'solid', 'blank')
LOSS$color = as.factor(as.numeric(as.factor(LOSS$interaction)))
LOSS$color = ifelse(LOSS$model == 'COAT', NA, LOSS$color)
g <- LOSS %>% subset(model!='COAT' ) %>%
  mutate(Var1 = as.character(Var1)) %>%
  mutate(Var1 = sapply(Var1,function(x){if(x == 'L_inf'){return('L_infinity')} else{return(as.character(x))}}),
         lambda = as.factor(lambda),
         mu_var = as.factor(mu_var),
         Model = paste0(model,', lambda = ',as.factor(lambda))) %>%
  #mutate(model = ifelse(model == 'COAT',model,paste0(model,', lambda = ',lambda))) %>%
  ggplot(aes(x = mu_var, y = value, color = Model,linetype = linetype)) + 
  geom_point(size=0.7) +
  scale_color_manual(values = c('red','blue','green','orange')) +
  facet_wrap(~Var1~Var2, scales = 'free_y', ncol = 4) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) +
  scale_linetype_identity() +
  geom_hline(data = loss, aes(yintercept = value,color='COAT')) +
  theme_bw() +
  labs(x = 'c', y = 'avg. loss') 
g <- g + guides(colour = guide_legend(override.aes = list( linetype = c(1,NA,NA,NA), shape = c(NA,16,16,16))))
ggsave(plot = g, width = 7, height = 7, dpi = 300, filename = paste0(result_dir,'clrcorr_loss.png'))


