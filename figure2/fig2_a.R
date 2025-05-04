library(tidyverse)
library(gdata)
library(ggpubr)
library(plyr)
library(ggsignif)
###plot
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
load("data/human_matched_out_plain.rda")
load("data/mouse_matched_out_plain.rda")
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat,diag = TRUE)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

human=flattenCorrMatrix(human_matched_out)
colnames(human)=c("sample1","sample2","human_cor")

mouse=flattenCorrMatrix(mouse_matched_out)
colnames(mouse)=c("sample1","sample2","mouse_cor")
infor=read.csv("../figure1/data//ribobase_7_17_23_manuscript_version.csv")
temp_plot <- function(human,species="human"){
  tmp1=merge(human,infor,by.x = "sample1",by.y="experiment_alias")
  tmp2=merge(tmp1,infor,by.x = "sample2",by.y="experiment_alias")
  dif_cell_dif_study=tmp2 %>% filter(cell_line.x!=cell_line.y & study_id.x!=study_id.y)
  same_cell_dif_study=tmp2 %>% filter(cell_line.x==cell_line.y & study_id.x!=study_id.y)
  same_cell_same_study=tmp2 %>% filter(cell_line.x==cell_line.y & study_id.x==study_id.y)
  dif_cell_dif_study_cor=data.frame(cor=dif_cell_dif_study[,3],status="diff.cell-line diff. studies",species)
  same_cell_dif_study_cor=data.frame(cor=same_cell_dif_study[,3],status="same cell-line diff. studies",species)
  #same_cell_same_study_cor = data.frame(cor=same_cell_same_study$human_cor, status="same cell-line same study")
  df_human_ribo=rbind(dif_cell_dif_study_cor,same_cell_dif_study_cor)
  return(df_human_ribo)
}
temp_plot2<- function(human,species="human"){
  tmp1=merge(human,infor,by.x = "sample1",by.y="experiment_alias")
  tmp2=merge(tmp1,infor,by.x = "sample2",by.y="experiment_alias")
  same_samaple_same_study=tmp2 %>% filter(sample1==sample2 & study_id.x==study_id.y)
  dif_samaples_dif_study=tmp2 %>% filter(sample1!=sample2 & study_id.x!=study_id.y)
  dif_samaples_same_study=tmp2 %>% filter(sample1!=sample2 & study_id.x==study_id.y)
  same_samaple_same_study_cor=data.frame(cor=same_samaple_same_study[,3],status="same sample same study",species)
  dif_samaples_dif_study_cor=data.frame(cor=dif_samaples_dif_study[,3],status="diff. samples diff. studies",species)
  dif_samaples_same_study_cor = data.frame(cor=dif_samaples_same_study[,3], status="diff. samples same study",species)
  df_human_ribo=rbind(same_samaple_same_study_cor,dif_samaples_dif_study_cor,dif_samaples_same_study_cor)
  return(df_human_ribo)
}

human_df <- temp_plot2(human,species="human")
mouse_df <- temp_plot2(mouse,species="mouse")
#pairwise.wilcox.test(human_df$cor, human_df$status, p.adjust.method = "BH")
#pairwise.wilcox.test(mouse_df$cor, mouse_df$status, p.adjust.method = "BH")
ribo_df <- rbind(human_df,mouse_df)
colnames(ribo_df) <- c("cor","status","Type")
ribo_df$Type=sub("human","Human",ribo_df$Type)
ribo_df$Type=sub("mouse","Mouse",ribo_df$Type)
ribo_df_fin <- ribo_df %>% unite("fin_group", Type:status, remove = FALSE,sep=":")
ribo_df_summary <- ribo_df_fin %>% dplyr::group_by(status,Type) %>% dplyr::summarise(round(mean(cor), digits=2))
colnames(ribo_df_summary) <- c("status","Type","cor_mean")
label_size = 10
legedn_size = 10

xx <- ddply(data.frame(ribo_df),.(Type,status),summarize,
            ymin = min(cor),
            ymax = max(cor),
            middle = median(cor),
            lower = quantile(cor,0.25),
            upper = quantile(cor,0.75))

ribo_df_plot <- ribo_df %>% mutate(anno=paste(status,Type))
#means <- aggregate(cor ~  status, ribo_df, FUN = function(x) {round(mean(x), digits = 3)})
pd = position_dodge(width = 0.5)
ggplot(ribo_df, aes(y=cor, x=Type,fill=status)) + 
  geom_boxplot(aes(alpha=status),width=0.3,outlier.shape = NA,position=pd,linetype="dotted",coef=1.5) + 
  stat_boxplot(width=0.3,aes(alpha=status,ymin = after_stat(lower), ymax = after_stat(upper)), outlier.shape = NA,position=pd,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = after_stat(ymax),alpha=status),position=pd,width=0.2,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = after_stat(ymin),alpha=status),position=pd,width=0.2,coef=1.5) +
   theme_bw() + 
  scale_fill_manual(values=c("#D8A0A8","#70A3C4","#DAA520"))+
  scale_alpha_manual(values = c(1,1,1))+
  guides(fill=guide_legend(nrow=1, byrow=TRUE),alpha=guide_legend(nrow=2, byrow=TRUE)) +
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
    legend.text = element_text(size=label_size))+labs(y= "R2 between RIBO and RNA",x="")
ggsave("ribo_rna_paired_r2.pdf", width = 2.4, height = 2.2, units = "in")

ggsave("ribo_rna_paired_r2_anno.pdf", width = 6.4, height = 2.2, units = "in")


