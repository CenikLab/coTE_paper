###paxDB_cor
library(randomcoloR)
library(edgeR)
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
#sample_infor=read.csv("../figure1/data/ribobase_7_17_23_manuscript_version.csv")
infor=read.csv("data/paxdb-uniprot-links-v4.2.tsv",header=T,sep="\t")
pax_cl = c("HepG2","A549","HEK293","HeLa","K562","MCF7","U2OS")

#out_fin <- data.frame(NULL)
df_cor <- function(pax_cl,df,type){
  out <- data.frame(matrix(nrow=7, ncol=1))
  rownames(out)=pax_cl
  for(i in pax_cl){
    print(i)
    #sample_infor_tmp <- sample_infor %>% filter(cell_line==i) %>%
     # dplyr::select(experiment_alias,cell_line,study_id)
    cl=read.csv(paste("data/9606-iBAQ_",i,"_Geiger_2012_uniprot.txt",sep=""),header=T,sep="\t")
    cl_anno=merge(cl,infor,by.x="string_external_id",by.y="ID")
    cl_anno[,4] <- sub("_HUMAN", "", cl_anno[,4])
    cl_df_fin=merge(cl_anno,df,by.x="gene_name",by.y=0)
    cl_df_fin=na.omit(cl_df_fin)
    print(dim(cl_df_fin))
    m <- cor(cl_df_fin[,i], cl_df_fin[,"abundance"],method="spearman")  
    out[i,1] <-as.numeric(m)  
  }
  #colnames(out)=pax_cl
  #write.csv(out,paste(i,type,"paxdb_cor_spearman.csv",sep="_"))
  return(out)
}

##flatten ribo
TE_flatten <- read.csv("data/human_TE_cellline_all_plain.csv",row.names=1)
TE = t(TE_flatten)
TE_prot_f <- TE[,pax_cl]
flatten_ribo_cor=df_cor(pax_cl,TE_prot_f,"TE_flatten")
colnames(flatten_ribo_cor)="cor"
flatten_ribo_cor$type="flatten"
flatten_ribo_cor$name=rownames(flatten_ribo_cor)

label_size = 8
legedn_size = 8
colors <- c("#D32F2F", "#1976D2", "#388E3C", "#FBC02D", "#8E44AD", "#E67E22", "#2C3E50")

ggplot(flatten_ribo_cor, aes(x = type, y = cor)) + 
  geom_boxplot(width = 0.2, outlier.shape = NA, linetype = "dotted", coef = 1.5) + 
  stat_boxplot(width = 0.2, aes(ymin = after_stat(lower), ymax = after_stat(upper)), outlier.shape = NA, coef = 1.5) +
  stat_boxplot(geom = "errorbar", aes(ymin = after_stat(ymax)), width = 0.1, coef = 1.5) +
  stat_boxplot(geom = "errorbar", aes(ymax = after_stat(ymin)), width = 0.1, coef = 1.5) +
  geom_point(aes(group = name, color = name), size = 1.5) +
  geom_line(aes(group = name, color = name)) +
  geom_text(aes(label = name, color = name), hjust = -0.1, vjust = 0.3, size = 2.5, show.legend = FALSE) +  # <-- this line adds labels
  theme_bw() + ylim(0.3,0.55) + 
  scale_color_manual(values = colors) +
  guides(fill = FALSE) +
  theme(
    axis.text = element_text(size = label_size),
    axis.title = element_text(size = label_size),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none",  # Remove legend since names are now on plot
  ) +
  labs(y = "Spearman correlation", x = "")

ggsave("TE_proteomics.pdf", width = 2, height = 2.2, units = "in")
