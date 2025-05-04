library(ComplexHeatmap)
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure4/")
df=read.csv("data/all_select_gene_global_new.csv",row.names=1)
load("../figure3/data/human_TE_rho.rda")
load("data/human_RNA_rho_new.rda")
#rownames(human_TE_rho) = colnames(human_TE_rho) = gsub("-", ".", rownames(human_TE_rho))
df_tmp <- df %>% filter(AUROC_TE > 0.8 & dff > 0.1)
#write.csv(df_tmp,"df_human_AUROC.csv")
df_map <- df_tmp %>% arrange(-AUROC_TE)
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat,diag = TRUE)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}
mapfunc <- function(i,human_TE_rho){
  l_tmp=strsplit(df_map[i,"gene_list"],",")[[1]]
  l = gsub("-", ".", l_tmp)
  print(l)
  tmp <- human_TE_rho[l,l]
  cor <- flattenCorrMatrix(tmp) %>% filter(cor!=1) %>% 
    mutate(abscor=abs(cor)) %>% arrange(-abscor)
  #tmp2 <- cor$cor
  #if (length(tmp2) < 231) {
    #tmp2<- c(tmp2, rep(0, 231 - length(tmp2)))
  #}
  return(cor$abscor)
}
lt = lapply(1:nrow(df_map), mapfunc,human_TE_rho=human_TE_rho)
lt_RNA = lapply(1:nrow(df_map), mapfunc,human_TE_rho=human_RNA_rho_new)

lt_new <- vector("list", length(lt) + length(lt_RNA))
lt_new[c(TRUE, FALSE)] <- lt
lt_new[c(FALSE, TRUE)] <- lt_RNA

df_col <- df_tmp %>% arrange(-dff) %>%
  dplyr::select(AUROC_TE,AUROC_RNA,anno,number_genes) %>% 
  pivot_longer(-c(anno,number_genes),names_to = "Type",values_to = "AUROC")
col_fun=paletteer::paletteer_d("rcartocolor::Temps")
library(circlize)
col_fun2=colorRamp2(c(10, 15, 20,25), c("#69D2E7FF","#F4B95AFF","#F38630FF","#FA6900FF"))

df_copy <- df_col
df_copy$anno[duplicated(df_copy$anno)] <- ""

ht_list = HeatmapAnnotation(AUROC = anno_barplot(df_col$AUROC, height = unit(1, "cm"),
                                                 gp = gpar(border=NA,fill = c("#984ea3","grey"),lty="blank"),
                                                axis_param = list(at = c(0, 0.2, 0.4, 0.6, 0.8), 
                                                                  labels = c("0%","20%", "40%", "60%", "80%"))),
                            annotation_name_rot = 90,gp = gpar(fontsize = 6)) %v%
  densityHeatmap(lt_new,show_quantiles = FALSE,title=NA,ylab="",
                 height = unit(1.5, "cm"),col = col_fun,range=c(0,0.6),
                 tick_label_gp = gpar(fontsize = 8)) %v%
  HeatmapAnnotation(size = df_col$number_genes,annotation_name_rot = 90,
                    gp = gpar(fontsize = 8,col = "black"),
                    simple_anno_size = unit(0.25, "cm"),
                    col=list(size=col_fun2)) %v%
  HeatmapAnnotation(foo=anno_text(df_copy$anno,
                                  gp = gpar(fontsize = 6),
                                  location = 1, rot = 55, 
                                  just = "right"))

draw(ht_list)
pdf("top_diff_update_4_12_24.pdf",width = 8,height = 3.5)
draw(ht_list)
dev.off()

