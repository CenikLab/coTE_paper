setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
df=read.csv("data/top5_data_update.csv")
label_size=8
df$rank=fct_rev(factor(df$rank, levels = df$rank))
ggplot(df,aes(y=rank, x=LOD))+ 
  geom_segment( aes(yend=rank, xend=0)) +
  geom_point(aes(color=species), size=3) +
  scale_y_discrete(labels=c("1" = "extracellular matrix structural constituent",
                            "2" = "extracellular matrix",
                            "3" = "collagen-containing extracellular matrix",
                            "4"="RNA binding",
                            "5"="extracellular space"))+
  theme_bw() + scale_color_manual(values=c("#984ea3","#4daf4a"))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "",x="LOD")
ggsave("mead_dif_GO_3_28_24.pdf",width = 4.5, height = 2.2, units = "in")
