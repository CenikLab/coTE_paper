setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure4/")

human_random_TE <- read.csv("data/auroc_GO_human_TE_random.csv",row.names=1)[,1,drop=FALSE]
human_random_RNA <- read.csv("data/auroc_GO_human_RNA_random.csv",row.names=1)[,1,drop=FALSE]
human_TE <- read.csv("data/auroc_te_GO_12_22_human.csv",row.names=1)[,1,drop=FALSE]
human_RNA <- read.csv("data/auroc_rna_GO_cl_human.csv",row.names=1)[,1,drop=FALSE]
tmp1 <- merge(human_random_TE,human_random_RNA,by=0)
colnames(tmp1) <- c("gene","random_TE_auroc","random_RNA_auroc")
tmp2 <- merge(human_TE,human_RNA,by=0)
colnames(tmp2) <- c("gene","TE_auroc","RNA_auroc")
humap_plot <- merge(tmp1,tmp2,by="gene") %>% pivot_longer(-gene,names_to = "ori_col",
                                                          values_to = "auroc")
humap_plot_fin <- humap_plot %>% mutate(type=if_else(ori_col=="random_TE_auroc" |
                                                ori_col=="TE_auroc","TE","RNA")) %>%
  mutate(status=if_else(ori_col=="random_TE_auroc" |
                   ori_col=="random_RNA_auroc","random","database")) %>%
  pivot_wider(-ori_col,names_from = type,values_from = auroc) %>% 
  mutate(species="human") %>% 
  mutate(type=paste(status,species))

"#CEF2E8","#F2DC9B"
label_size=8
ggplot(humap_plot_fin,aes(x=RNA,y=TE,color=type))+
  geom_point(size=0.5,alpha=0.4)+
  theme_bw() +scale_color_manual(values=c("#4174D9","#F28705"))+
  geom_rug(data=humap_plot_fin[humap_plot_fin$status=="random",],
           aes(x=RNA,y=TE),alpha=0.3,sides = "tr")+
  geom_rug(data=humap_plot_fin[humap_plot_fin$status=="database",],
           aes(x=RNA,y=TE),alpha=0.3,,sides = "tr",outside=TRUE)+
  coord_cartesian(clip = "off")+
  geom_hline(yintercept = 0.7, linetype="dashed", 
             color = "grey", size=1)+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size))+labs(y= "TE AUROC",x="RNA AUROC")
ggsave("GO_cor_random_human.pdf", width = 2.4, height = 2.2, units = "in")
