library(tidyverse)
label_size = 10
legedn_size = 10
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure1/")
human_df=read.csv("data/human_QC_summary.csv")
human_df_tmp <- human_df %>% mutate(species="human") %>% mutate(type=if_else(coverage<0.1,"<0.1X",
                                                                             if_else(coverage<1,"<2X",
                                                                                     if_else(coverage<5,"<5X",">=5X")))) %>%
  select(-status)

mouse_df=read.csv("data/mouse_QC_summary_fin.csv")
mouse_df_tmp <- mouse_df %>% mutate(species="mouse") %>% mutate(type=if_else(coverage<0.1,"<0.1X",
                                                                             if_else(coverage<1,"<2X",
                                                                                     if_else(coverage<5,"<5X",">=5X"))))
df_all <- rbind(human_df_tmp,mouse_df_tmp)

df_group_summary <- df_all %>% group_by(species,type) %>% summarise(count=n()) 
label_size = 10
legedn_size = 10
ggplot(df_group_summary,aes(y=type,x=count,fill=species)) +geom_col(position = "dodge",width=0.8) +
  theme_bw() +
  scale_fill_manual(values=c("#984ea3","#4daf4a"))+
  guides(fill=guide_legend(nrow=1, byrow=TRUE),alpha=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_discrete(labels=c("<0.1X" = "0-0.1", "<2X" = "0.1-2", 
                            "<5X" = "2-5", ">=5" = ">=5"))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size))+labs(y= "Coverage (cutoff: 0.1X)",
                                                          x="")
ggsave("riboseq_Yue_QC_Coverage_0.1_fin.pdf", width = 2.4, height = 2.2, units = "in")

