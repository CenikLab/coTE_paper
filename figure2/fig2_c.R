library(propr)
library(tidyverse)
library(gdata)
library(ggpubr)

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}


setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
load("data/human_ribo_cor_flatten.rda")
load("data/mouse_ribo_cor_flatten.rda")
###ribo
human_ribo=flattenCorrMatrix(human_ribo_cor)
colnames(human_ribo)=c("sample1","sample2","human_ribo_cor")

mouse_ribo=flattenCorrMatrix(mouse_ribo_cor)
colnames(mouse_ribo)=c("sample1","sample2","mouse_ribo_cor")

###RNA
load("data/human_rna_cor_flatten.rda")
load("data/mouse_rna_cor_flatten.rda")
###ribo
human_rna=flattenCorrMatrix(human_rna_cor)
colnames(human_rna)=c("sample1","sample2","human_rna_cor")

mouse_rna=flattenCorrMatrix(mouse_rna_cor)
colnames(mouse_rna)=c("sample1","sample2","mouse_rna_cor")


#infor=read.csv("infor.csv")
infor_tmp=read.csv("../figure1/data//ribobase_7_17_23_manuscript_version.csv")
infor = infor_tmp %>%  dplyr::select(experiment_alias,cell_line,study_id,organism) 
colnames(infor) = c("experiment_alias","cell_line","study_id","species")
temp_plot <- function(human_ribo,species="human"){
  tmp1=merge(human_ribo,infor,by.x = "sample1",by.y="experiment_alias")
  tmp2=merge(tmp1,infor,by.x = "sample2",by.y="experiment_alias")
  dif_cell_dif_study=tmp2 %>% filter(cell_line.x!=cell_line.y & study_id.x!=study_id.y)
  same_cell_dif_study=tmp2 %>% filter(cell_line.x==cell_line.y & study_id.x!=study_id.y)
  same_cell_same_study=tmp2 %>% filter(cell_line.x==cell_line.y & study_id.x==study_id.y)
  dif_cell_dif_study_cor=data.frame(cor=dif_cell_dif_study[,3],status="diff.cell-line diff. studies",species)
  same_cell_dif_study_cor=data.frame(cor=same_cell_dif_study[,3],status="same cell-line diff. studies",species)
  #same_cell_same_study_cor = data.frame(cor=same_cell_same_study$human_ribo_cor, status="same cell-line same study")
  df_human_ribo=rbind(dif_cell_dif_study_cor,same_cell_dif_study_cor)
  return(df_human_ribo)
}
human_ribo_df <- temp_plot(human_ribo,species="Human_ribo")
mouse_ribo_df <- temp_plot(mouse_ribo,species="Mouse_ribo")

human_rna_df <- temp_plot(human_rna,species="Human_rna")
mouse_rna_df <- temp_plot(mouse_rna,species="Mouse_rna")
ribo_df <- rbind(human_ribo_df,mouse_ribo_df,human_rna_df,mouse_rna_df)
ribo_df_summary <- ribo_df %>% dplyr::group_by(status,species) %>% dplyr::summarise(med=round(median(cor), digits=2))
wilcox.test(cor ~ status, data = human_ribo_df)
wilcox.test(cor ~ status, data = human_rna_df)
wilcox.test(cor ~ status, data = mouse_ribo_df)
wilcox.test(cor ~ status, data = mouse_rna_df)
label_size = 8
legedn_size = 8
#means <- aggregate(cor ~  status, ribo_df, FUN = function(x) {round(mean(x), digits = 3)})
pd = position_dodge(width = 0.5)
ggplot(ribo_df, aes(y=cor, x=species,fill=status)) + 
  geom_boxplot(aes(alpha=status),width=0.2,outlier.shape = NA,position=pd,linetype="dotted",coef=1.5) + 
  stat_boxplot(width=0.2,aes(alpha=status,ymin = after_stat(lower), ymax = after_stat(upper)), outlier.shape = NA,position=pd,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = after_stat(ymax),alpha=status),position=pd,width=0.2,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = after_stat(ymin),alpha=status),position=pd,width=0.2,coef=1.5) +
  theme_bw() + 
  scale_fill_manual(values=c("#D8A0A8","#70A3C4"))+
  scale_alpha_manual(values = c(1,1))+
  guides(fill=guide_legend(nrow=2, byrow=TRUE),alpha=guide_legend(nrow=2, byrow=TRUE)) +
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "Spearman correlation",x="")
ggsave("r2_spearman.pdf", width = 5, height = 2.2, units = "in")
