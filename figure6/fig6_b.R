library(tidyverse)
library(plyr)
library(ggpubr)

label_size = 8
legedn_size = 8
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure6/")
human_data=read.csv("data/auroc_te_humap.csv")
colnames(human_data)=c("term","mean","median","uq","number")
human_data[,'Type']='TE'
human_data[,'Type2']='RiboBase'
anno=read.csv("data/humap_TE_RNA_gene.csv",row.names=1)[,c(1,5)]
df_anno <- anno %>% distinct(gene_list, .keep_all = TRUE)
human_data_fin <- human_data %>% filter(term %in% df_anno$term)

human_bg=read.csv("data/auroc_humap_TE_random.csv")
colnames(human_bg)=c("term","mean","median","uq","number")
human_bg[,'Type']='TE'
human_bg[,'Type2']='Random'
human_bg_fin <- human_bg %>% filter(term %in% df_anno$term)

human_data_rna=read.csv("data/auroc_rna_humap_cl.csv")
colnames(human_data_rna)=c("term","mean","median","uq","number")
human_data_rna[,'Type']='RNA'
human_data_rna[,'Type2']='RiboBase'
human_data_rna_fin <- human_data_rna %>% filter(term %in% df_anno$term)

human_bg_rna=read.csv("data/auroc_humap_RNA_random_cl.csv")
colnames(human_bg_rna)=c("term","mean","median","uq","number")
human_bg_rna[,'Type']='RNA'
human_bg_rna[,'Type2']='Random'
human_bg_rna_fin <- human_bg_rna %>% filter(term %in% df_anno$term)

tmp=merge(human_data_fin,human_data_rna,by="term")
#write.csv(tmp,"TE_RNA_positive.csv")

all <- rbind(human_data_fin,human_data_rna_fin,human_bg_fin,human_bg_rna_fin)
all_plot <- all %>% dplyr::select(term,mean,number,Type,Type2)

all_plot$type_fin <- paste(all_plot$Type, all_plot$Type2, sep = "_")
#write.csv(all_plot,"humap_AUROC_sig.csv")
#all_plot=read.csv("humap_AUROC_sig.csv")
combination=combn(unique(all_plot$type_fin), 2)
for (i in 1:ncol(combination)){
  tmp=combination[,i]
  wilcox=wilcox.test(all_plot[all_plot$type_fin==tmp[1],]$mean,all_plot[all_plot$type_fin==tmp[2],]$mean,paired=TRUE)
  wilcox_p=wilcox$p.value
  print(tmp)
  print(wilcox_p)
}

fin_summary <- all_plot %>% dplyr::group_by(Type,Type2) %>% 
  dplyr::summarise(auroc=round(median(mean), digits=2))

###median
pd = position_dodge(width = 0.8)
ggplot(all_plot, aes(x=mean)) + 
  geom_histogram(data=all_plot[all_plot$Type=="TE",],aes(x=mean,fill=Type2),
                 position = "identity", bins = 30,alpha=0.7)+
  geom_histogram(data=all_plot[all_plot$Type=="RNA",],
                 aes(x=mean,color=Type2),fill="transparent",
                 position = "identity", bins = 30)+
  theme_bw() + 
  scale_fill_manual(values=c("#CEF2E8","#F2DC9B"))+
  scale_color_manual(values=c("#4174D9","#F28705"))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "AUROC",x="")
ggsave("humap_AUROC_hist.pdf", width = 2.4, height = 2.2, units = "in")

p1=all_plot %>% filter(Type=="RNA"&Type2=="RiboBase")
by <- 0.03
h1 <- ggplot(data = p1, aes(x = mean)) +
  geom_histogram(breaks = seq(from = 0, to = 1, by = by),
                 color = "red", fill = "transparent")
df1 <- ggplot_build(h1)$data[[1]][ , c("x", "y")]

p2=all_plot %>% filter(Type=="RNA"&Type2=="Random")
h2 <- ggplot(data = p2, aes(x = mean)) +
  geom_histogram(breaks = seq(from = 0, to = 1, by = by),
                 color = "red", fill = "transparent")
df2 <- ggplot_build(h2)$data[[1]][ , c("x", "y")]

ggplot() + 
  geom_histogram(data=all_plot[all_plot$Type=="TE",],aes(x=mean,fill=Type2),
                 position = "identity", bins = 30,alpha=0.75,
                 breaks = seq(from = 0, to = 1, by = by))+
  geom_step(data = df1, aes(x = x - by/2, y = y),color="#F28705")+
  geom_step(data = df2, aes(x = x - by/2, y = y),color="#4174D9")+
  theme_bw() + 
  scale_fill_manual(values=c("#CEF2E8","#F2DC9B"))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "",x="AUROC")
ggsave("humap_AUROC_hist.pdf", width = 4, height = 2.2, units = "in")
