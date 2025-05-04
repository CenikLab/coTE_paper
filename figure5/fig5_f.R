library(ggrepel)
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure5/")
###df for figure f
df=read.csv("data/LRRC_multimer2.csv")
label_size=8
df_28=df %>% filter(Type=="LRRC28")
ggplot(df_28,aes(x=iptm.tpm,y=pDOCKQ)) + geom_point(color="#BA93DF")+
  geom_text_repel(aes(label=gene), size=2,
                  box.padding = 0.5, 
                  point.padding = 0.5, 
                  segment.color = 'grey50',max.overlaps = Inf) +
  theme_bw() + geom_vline(xintercept = 0.7, linetype="dotted", 
                          color = "grey", size=0.5)+
  geom_hline(yintercept = 0.23, linetype="dotted", 
             color = "grey", size=0.5)+ xlim(0.1, 0.8)+ylim(0,0.5)+
  theme(axis.text=element_text(size=label_size,colour = "black"),
        axis.title = element_text(size=label_size,colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="left",
        legend.text = element_text(size=label_size,face="bold")) +
  labs(x= "ipTM+pTM",y="pDOCKQ")
ggsave("alphafold_TF_28.pdf",width =2.5 , height = 2.2, units = "in")

##compare with LRRC42
df_42=df %>% filter(Type=="LRRC42")
ggplot(df_42,aes(x=iptm.tpm,y=pDOCKQ)) + geom_point(color="grey")+
  geom_text_repel(aes(label=gene), size=2,
                  box.padding = 0.5, 
                  point.padding = 0.5, 
                  segment.color = 'grey50',max.overlaps = Inf) +
  theme_bw() + geom_vline(xintercept = 0.7, linetype="dotted", 
                          color = "grey", size=0.5)+
  geom_hline(yintercept = 0.23, linetype="dotted", 
             color = "grey", size=0.5)+
  xlim(0.1, 0.8)+ylim(0,0.5)+
  theme(axis.text=element_text(size=label_size,colour = "black"),
        axis.title = element_text(size=label_size,colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="left",
        legend.text = element_text(size=label_size,face="bold")) +
  labs(x= "ipTM+pTM",y="pDOCKQ")
ggsave("alphafold_TF_42.pdf",width =2.5 , height = 2.2, units = "in")

###df for figure e
df=read.csv("data/LRRC_TF_value.csv")
####scatter plot FOX family
outlier_threshold_x <- 0.7
outlier_threshold_y <- 0.23
df$is_outlier <- ifelse(df$iptm.tmp > outlier_threshold_x | df$pDOCKQ > outlier_threshold_y, TRUE, FALSE)
ggplot(df,aes(x=iptm.tmp,y=pDOCKQ)) + geom_point(color="darkgreen",size=0.8)+
  geom_text_repel(data = subset(df, is_outlier), aes(label = gene), size=2,
                  box.padding = 0.5, 
                  point.padding = 0.5, 
                  segment.color = 'grey50',max.overlaps = Inf) +
  theme_bw() + geom_vline(xintercept = 0.7, linetype="dotted", 
                          color = "grey", size=0.5)+
  geom_hline(yintercept = 0.23, linetype="dotted", 
             color = "grey", size=0.5)+ xlim(0.2,0.8)+ylim(0.05,0.5)+
  theme(axis.text=element_text(size=label_size,colour = "black"),
        axis.title = element_text(size=label_size,colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="left",
        legend.text = element_text(size=label_size,face="bold")) +
  labs(x= "ipTM+pTM",y="pDOCKQ")
ggsave("alphafold_FOX_family.pdf",width =2.5 , height = 2.2, units = "in")



