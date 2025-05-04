setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
homo_gene <- read.csv("data/homo_gene_ian.csv")

human_TE <- read.csv("data/human_TE_cellline_all_plain.csv",row.names=1)
mouse_TE <- read.csv("data/mouse_TE_cellline_all_plain.csv",row.names=1)

#colnames(human_TE) <- tolower(colnames(human_TE))
#colnames(mouse_TE) <- tolower(colnames(mouse_TE))

human_TE_homo <- human_TE[,homo_gene$human]
mouse_TE_homo <- mouse_TE[,homo_gene$mouse]
colnames(human_TE_homo) <- homo_gene$marker
colnames(mouse_TE_homo) <- homo_gene$marker
label_size = 10
legedn_size = 10

human_TE_homo_tmp <- t(human_TE_homo)
human_TE_homo_mean <- data.frame(mean_TE=apply(human_TE_homo_tmp, 1, mean))

mouse_TE_homo_tmp <- t(mouse_TE_homo)
mouse_TE_homo_mean <- data.frame(mean_TE=apply(mouse_TE_homo_tmp, 1, mean))

fin_plot <- merge(human_TE_homo_mean,mouse_TE_homo_mean,by=0)
colnames(fin_plot) <- c("gene","human_TE_mean","mouse_TE_mean")
cor(fin_plot$human_TE_mean,fin_plot$mouse_TE_mean,method="spearman")

model <- lm(human_TE_mean ~ mouse_TE_mean, data = fin_plot)
conf.int <- predict(model, interval = "prediction")
mydata <- cbind(fin_plot, conf.int)
#mydata_plot <- mydata %>%
  #mutate(highlight=if_else(human_TE_mean>upr|human_TE_mean<lwr,
#                                        "yes","no"))

#mydata=read.csv("data/fin_human_TE.csv")
mydata_plot <- mydata %>%
  mutate(highlight=if_else(human_TE_mean>upr,
                           "human",if_else(
                             human_TE_mean<lwr,"mouse","no")))

p_intervel <- ggplot(data = mydata_plot, aes(x = mouse_TE_mean, y = human_TE_mean, color=highlight)) +
  geom_point(size=0.1) +
  geom_line(aes(y = lwr), color = "#48736F", linetype = "dashed")+
  geom_line(aes(y = upr), color = "#48736F", linetype = "dashed") +
  scale_color_manual(values=c("#984ea3","#4daf4a","grey"))+
  annotate("text", x = Inf, y = Inf, label = "Spearman: 0.90", vjust = 1, hjust = 1, size = 3, color = "black")+
theme_bw() +scale_y_continuous(limits=c(-3, 3), breaks=seq(-3, 3, by=1))+
  scale_x_continuous(limits=c(-3, 3), breaks=seq(-3, 3, by=1))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "Human avg. TE",x="Mouse avg. TE")
mydata_detail <- mydata_plot %>%
  mutate(highlight=if_else(human_TE_mean>upr,
                           "human_high",if_else(human_TE_mean<lwr,
                                                "human_low","no")))
#write.csv(mydata_detail,"mydata_plot_95_predictCI_ian.csv")

#mydata_plot=read.csv("data/fin_human_TE.csv")
p <- ggplot(data = mydata_plot, aes(x = mouse_TE_mean, y = human_TE_mean)) +
  geom_point(size=0.1,col="#666362",alpha=0.8) +
  annotate("text", x = Inf, y = Inf, label = "Spearman: 0.90", vjust = 1, hjust = 1, size = 3, color = "black")+
  theme_bw() + 
  scale_y_continuous(limits=c(-3, 3), breaks=seq(-3, 3, by=1))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "Human avg. TE",x="Mouse avg. TE")
pdf("TE_cor_homo_ci_ian.pdf",
    width = 2.4, height = 2.2)
ggExtra::ggMarginal(p_intervel, type = "histogram",
                    xparams = list(fill = "#4daf4a"),
                    yparams = list(fill = "#984ea3"))
dev.off()
ggsave("TE_cor_homo.pdf", width = 2.4, height = 2.2, units = "in")

mydata_detail$human <- toupper(mydata_detail$gene)
mydata_detail$mouse <- str_to_title(mydata_detail$gene)
write.csv(mydata_detail,"mydata_plot_95_predictCI.csv")
