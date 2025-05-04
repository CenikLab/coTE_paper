library(tidyverse)
label_size = 10
legedn_size = 10
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure1/")
human_df=read.csv("data/human_QC_summary.csv")
human_df_tmp <- human_df %>% mutate(type="human") %>% select(-status)
mouse_df=read.csv("data/mouse_QC_summary_fin.csv")
mouse_df_tmp <- mouse_df %>% mutate(type="mouse")

df = rbind(human_df_tmp,mouse_df_tmp)

ggplot(df,aes(y=CDS_pct,x=type,fill=type))+geom_violin(alpha=0.5) +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  geom_hline(yintercept=0.7, linetype="dashed", color = "#C07A92", size=2) +
  theme_bw() +
  scale_fill_manual(values=c("#984ea3","#4daf4a"))+
  guides(fill=guide_legend(nrow=1, byrow=TRUE),alpha=guide_legend(nrow=2, byrow=TRUE)) +
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "CDS percentage (cutoff: 0.7)",
                                                                      x="")
ggsave("riboseq_Yue_QC_CDS_0.7_fin.pdf", width = 2.4, height = 2.2, units = "in")
