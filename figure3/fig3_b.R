setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure3/")
library(ggplot2)
label_size=6
label_size=8
###the bin files: how many gene pairs are located in the certain range. you can summaries yourself!
load("data/human_mouse_TE_bin.rda")
load("data/human_mouse_TE_r_bin.rda")
bin_size=0.1
# Create scatter plot with bins
cor(as.numeric(human_mouse_TE_bin$cor1_bin),as.numeric(human_mouse_TE_bin$cor2_bin),method="spearman")
cor(as.numeric(human_mouse_TE_r_bin$cor1_bin),as.numeric(human_mouse_TE_r_bin$cor2_bin),method="spearman")
ggplot(human_mouse_TE_bin, aes(x=as.numeric(as.character(cor1_bin)), y=as.numeric(as.character(cor2_bin)),
                              size=log10(human), color=log10(human))) + 
  geom_point(alpha=0.8) + 
  scale_size_continuous(range = c(0.1, 2)) +
  scale_color_gradient2(low="black", mid="lightblue", high="darkred", midpoint=log10(median(human_mouse_TE_bin$human)))+ # Colors can be adjusted
  theme_bw() + 
  theme(axis.text=element_text(size=label_size,colour = "black"),
        axis.title = element_text(size=label_size,colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="left",
        legend.text = element_text(size=label_size,face="bold",colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=1))+
  labs( x="human",y="mouse",size="log10(count)")+
  scale_x_continuous(breaks=seq(-1, 1, by=bin_size), limits=c(-1, 1)) +
  scale_y_continuous(breaks=seq(-1, 1, by=bin_size), limits=c(-1, 1)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed")
ggsave("CoTE_r_new.pdf",width = 3.5, height = 3, units = "in")
