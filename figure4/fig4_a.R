library(tidyverse)
library(plyr)
library(ggpubr)

label_size = 10
legedn_size = 10

setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure4/")
human_data=read.csv("data/auroc_te_GO_12_22_human.csv")
colnames(human_data)=c("term","mean","median","uq","number")
human_data[,'Type']='TE'
human_data[,'Type2']='Human'

human_data_rna=read.csv("data/auroc_rna_GO_cl_human.csv")
colnames(human_data_rna)=c("term","mean","median","uq","number")
human_data_rna[,'Type']='RNA'
human_data_rna[,'Type2']='Human'

mouse_data=read.csv("data/auroc_te_GO_12_22_mouse.csv")
colnames(mouse_data)=c("term","mean","median","uq","number")
mouse_data[,'Type']='TE'
mouse_data[,'Type2']='Mouse'

mouse_data_rna=read.csv("data/auroc_rna_GO_cl_mouse.csv")
colnames(mouse_data_rna)=c("term","mean","median","uq","number")
mouse_data_rna[,'Type']='RNA'
mouse_data_rna[,'Type2']='Mouse'

all <- rbind(human_data,human_data_rna,mouse_data,mouse_data_rna)
all_plot <- all %>% dplyr::select(term,mean,number,Type,Type2) %>% 
  mutate(group=paste(Type,Type2,sep="_"))

fin_summary <- all_plot %>% dplyr::group_by(Type,Type2) %>% 
  dplyr::summarise(auroc=round(median(mean), digits=2))

all_plot$group <- factor(all_plot$group, levels = c("RNA_Human", "TE_Human", "RNA_Mouse", "TE_Mouse"))
###mean ture auroc
pd = position_dodge(width = 0.8)

p<-ggplot(all_plot, aes(y=mean, x=group,fill=group)) + 
  geom_boxplot(width=0.3,outlier.shape = NA,position=pd,linetype="dotted",coef=1.5) + 
  stat_boxplot(width=0.3,aes(ymin = after_stat(lower), ymax = after_stat(upper)),outlier.shape = NA, position=pd,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = after_stat(ymax)),position=pd,width=0.2,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = after_stat(ymin)),position=pd,width=0.2,coef=1.5) +
  theme_bw() +scale_fill_manual(values=c("#b3a2c7","#984ea3","#a8ddb5","#4daf4a"))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size))+labs(y= "AUROC for GO",x="")

ggsave("GO_AUROC_mean.pdf", width = 2.4, height = 2.2, units = "in")
