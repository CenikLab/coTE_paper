setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
library(compositions)
library(ppcor)
homo <- read.csv("data/homo_gene_ian.csv",row.names=1)

human_TE <- read.csv("data/human_TE_cellline_all_plain.csv",row.names=1)
mouse_TE <- read.csv("data/mouse_TE_cellline_all_plain.csv",row.names=1)

human_TE_homo <- human_TE[,homo$human]
mouse_TE_homo <- mouse_TE[,homo$mouse]

colnames(human_TE_homo) <- homo$marker
colnames(mouse_TE_homo) <- homo$marker

label_size = 8
legedn_size = 8

human_TE_homo_tmp <- t(human_TE_homo)
human_TE_homo_se <- data.frame(mean_TE=apply(human_TE_homo_tmp, 1, msd))
human_TE_homo_mean <- data.frame(mean_TE=apply(human_TE_homo_tmp, 1, function(x) mean(x,robust=TRUE)))
mouse_TE_homo_tmp <- t(mouse_TE_homo)
mouse_TE_homo_se <- data.frame(mean_TE=apply(mouse_TE_homo_tmp, 1, msd))
mouse_TE_homo_mean <- data.frame(mean_TE=apply(mouse_TE_homo_tmp, 1, function(x) mean(x,robust=TRUE)))

fin_plot_tmp1 <- merge(human_TE_homo_se,mouse_TE_homo_se,by=0)
colnames(fin_plot_tmp1) <- c("gene","human_TE_msd","mouse_TE_msd")
fin_plot_tmp2 <- merge(fin_plot_tmp1,human_TE_homo_mean,by.x="gene",by.y=0)
colnames(fin_plot_tmp2)<- c("gene","human_TE_msd","mouse_TE_msd","human_TE")
fin_plot <- merge(fin_plot_tmp2,mouse_TE_homo_mean,by.x="gene",by.y=0)
colnames(fin_plot)<- c("gene","human_TE_msd","mouse_TE_msd","human_TE","mouse_TE")
pcor.test(fin_plot$human_TE_msd,fin_plot$mouse_TE_msd,
          fin_plot[,c("human_TE","mouse_TE")],method="spearman")
X_resid<-resid(lm(fin_plot$human_TE_msd~fin_plot$human_TE + fin_plot$mouse_TE,fin_plot))
Y_resid<-resid(lm(fin_plot$mouse_TE_msd~fin_plot$human_TE + fin_plot$mouse_TE,fin_plot))
cor(X_resid, fin_plot$human_TE, method="spearman")
cor(Y_resid, fin_plot$mouse_TE, method="spearman")
cor(X_resid,Y_resid,method="spearman")
ggplot(data = fin_plot, aes(x = human_TE, y = X_resid)) +
  geom_point(size=0.05,col="#666362",alpha=0.8) + 
  theme_bw() + 
  annotate("text", x = Inf, y = Inf, label = "Spearman: -0.025", vjust = 1, hjust = 1, size = 3, color = "black")+
  theme(axis.text=element_text(size=label_size,colour = "black"),
        axis.title = element_text(size=label_size,colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold",colour = "black"))+
  labs(y= "Adjusted msd of human TE",x="mean of human TE")
ggsave("Human_MAD_vs_mean.pdf",
    width = 2.4, height = 2.2)

ggplot(data = fin_plot, aes(x = mouse_TE, y = Y_resid)) +
  geom_point(size=0.05,col="#666362",alpha=0.8) + 
  theme_bw() + 
  annotate("text", x = Inf, y = Inf, label = "Spearman: -0.033", vjust = 1, hjust = 1, size = 3, color = "black")+
  theme(axis.text=element_text(size=label_size,colour = "black"),
        axis.title = element_text(size=label_size,colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold",colour = "black"))+
  labs(y= "Adjusted msd of mouse TE",x="mean of mouse TE")
ggsave("mouse_MAD_vs_mean.pdf",
    width = 2.4, height = 2.2)

mydata_plot=read.csv("data/MSD_data.csv")
p_anno<- ggplot(data = mydata_plot, aes(x = X_resid, y = Y_resid, color=highlight)) +
  geom_point(size=0.1) +
  geom_line(aes(y = lwr), color = "#48736F", linetype = "dashed")+
  geom_line(aes(y = upr), color = "#48736F", linetype = "dashed") +
  scale_color_manual(values=c("#984ea3","#4daf4a","grey"))+
  annotate("text", x = Inf, y = Inf, label = "Spearman: 0.73", vjust = 1, hjust = 1, size = 3, color = "black")+
  theme_bw() +
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "Adjusted human TE msd",x="Adjusted mouse TE msd")

write.csv(mydata_detail,"mydata_plot_95_predictCI_msd_ian_3-72-24.csv")


p <- ggplot(data = mydata_plot, aes(x = X_resid, y = Y_resid)) +
  geom_point(size=0.1,col="#666362",alpha=0.8) +
  annotate("text", x = Inf, y = Inf, label = "Spearman: 0.61", vjust = 1, hjust = 1, size = 3, color = "black")+
  theme_bw() + 
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "Human avg. TE MAD",x="Mouse avg. TE MAD")
pdf("TE_mad_homo_ian_3-27-24.pdf",
    width = 2.4, height = 2.2)
ggExtra::ggMarginal(p_anno, type = "histogram",
                    xparams = list(fill = "#4daf4a"),
                    yparams = list(fill = "#984ea3"))
dev.off()
ggsave("TE_cor_homo.pdf", width = 2.4, height = 2.2, units = "in")
not_human=read.csv("/Users/yliu5/Library/CloudStorage/Box-Box/Canlab/publication_2023_version/Figure/Figure3/GO_bg/msd_not_list.csv")
mouse_select <- mydata_plot %>% filter(highlight=="mouse") %>% 
  filter(mouse %in% not_human)

df <- data.frame(x = X_resid, y = Y_resid,
                 d = densCols(X_resid, Y_resid, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

ggplot(df) +
  geom_point(aes(x, y, col = d), size = 1) +
  scale_color_identity() +
  theme_bw()


ggplot(df) +
  geom_point(aes(x, y, col = d), size = 0.5) +
  scale_color_identity() +
  annotate("text", x = Inf, y = Inf, label = "Spearman: 0.63", vjust = 1, hjust = 1, size = 3, color = "black")+
  theme_bw() +
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size,face="bold"))+labs(y= "Adjusted human TE msd",x="Adjusted mouse TE msd")
ggsave("TE_msd.pdf", width = 2.4, height = 2.2, units = "in")
