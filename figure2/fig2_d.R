flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

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
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
load("data/human_te_cor_flatten.rda")
load("data/mouse_te_cor_flatten.rda")

#infor=read.csv("infor.csv")
infor_tmp=read.csv("../figure1/data/ribobase_7_17_23_manuscript_version.csv")
infor = infor_tmp %>%  dplyr::select(experiment_alias,cell_line,study_id,organism) 
colnames(infor) = c("experiment_alias","cell_line","study_id","species")

##rna
human_TE=flattenCorrMatrix(human_te_cor)
colnames(human_TE)=c("sample1","sample2","human_rna_cor")

mouse_te=flattenCorrMatrix(mouse_te_cor)
colnames(mouse_te)=c("sample1","sample2","mouse_rna_cor")

human_te_df <- temp_plot(human_TE,species="Human")
mouse_te_df <- temp_plot(mouse_te,species="Mouse")
wilcox.test(cor ~ status, data = human_te_df)
wilcox.test(cor ~ status, data = mouse_te_df)
t.test(human_te_df[human_te_df$status=="diff.cell-line diff. studies",]$cor,
       human_te_df[human_te_df$status=="same cell-line diff. studies",]$cor)
te_df <- rbind(human_te_df,mouse_te_df)
colnames(te_df) <- c("cor","status","Type")

te_df_summary <- te_df %>% dplyr::group_by(Type,status) %>% dplyr::summarise(round(median(cor), digits=2))
colnames(te_df_summary) <- c("Type","status","cor_median")
label_size = 8
legedn_size = 8
#means <- aggregate(cor ~  status, ribo_df, FUN = function(x) {round(mean(x), digits = 3)})
pd = position_dodge(width = 0.5)
ggplot(te_df, aes(y=cor, x=Type,fill=status)) + 
  geom_boxplot(width=0.2,outlier.shape = NA,position=pd,linetype="dotted",coef=1.5) + 
  stat_boxplot(width=0.2,aes(ymin = after_stat(lower), ymax = after_stat(upper)), outlier.shape = NA,position=pd,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = after_stat(ymax)),position=pd,width=0.2,coef=1.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = after_stat(ymin)),position=pd,width=0.2,coef=1.5) +
  theme_bw() + coord_cartesian(ylim=c(0,1))+
  scale_fill_manual(values=c("#D8A0A8","#70A3C4"))+
  geom_tile(aes(y=NA_integer_, alpha = factor(status))) + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE),alpha=guide_legend(nrow=2, byrow=TRUE)) +
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "Spearman correlation",x="")
ggsave("TE_r2.pdf", width = 2.4, height = 2.2, units = "in")

