library(circlize)
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure6/")
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat,diag = TRUE)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

load("data/cor_TE_exocyst.rda")
cor_TE_exocyst[cor_TE_exocyst<0]=0
TE_d_tmp <- flattenCorrMatrix(cor_TE_exocyst)
TE_d <- TE_d_tmp %>% filter(cor!=1)

load("data/cor_RNA_exocyst.rda")
cor_RNA_exocyst[cor_RNA_exocyst<0]=0
RNA_d_tmp <- flattenCorrMatrix(cor_RNA_exocyst)
RNA_d <- RNA_d_tmp %>% filter(cor!=1)

tmp_clean <- merge(TE_d,RNA_d,by.x=c("row","column"),by.y=c("row","column"))
colnames(tmp_clean) <- c("row","column","TE_cor","RNA_cor")

net_edges <- tmp_clean %>% dplyr::select(row,column)
net_all_leaves <- colnames(cor_TE_exocyst)
clean_cor_tmp <- tmp_clean %>% filter(abs(TE_cor)>=0.1 | abs(RNA_cor)>=0.1)
clean_cor <- clean_cor_tmp %>% mutate(group=if_else(abs(TE_cor)>=0.1 & abs(RNA_cor)>=0.1,
                                                     "Both",if_else(
                                                       abs(TE_cor)>=0.1 &  abs(RNA_cor)<0.1,
                                                       "TE_only","RNA_only"
                                                     )))
clean_cor_fin <- clean_cor %>% 
  mutate(cor_color=if_else(group=="Both","#B9D5E9",
                         if_else(group=="TE_only","#984ea3","#CCCCCC")))


clean_cor_fin <- clean_cor_fin %>% 
  mutate(fin_cor=if_else(group=="Both",abs(TE_cor),
                           if_else(group=="TE_only",abs(TE_cor),abs(RNA_cor))))

fin_mat <- clean_cor_fin %>% arrange(cor_color)

input_mat <- fin_mat %>% dplyr::select("row","column","fin_cor")
pdf("exocyst_circle_cor.pdf",
    width =2.4, height = 2.4)
chordDiagram(input_mat,grid.col="grey",col=fin_mat$cor_color,
             transparency = 0.1,
             annotationTrack = "grid")
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(median(xlim), median(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE,cex=0.5)
}
dev.off()
