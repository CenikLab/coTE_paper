setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure3/")
gene=c("RPS3","RPS15A","RPS3A","HK2")

TE_df = read.csv("../figure2/data/human_TE_cellline_all_plain.csv",row.names=1)

TE_df_cotrans <- TE_df[,gene]

#ha = rowAnnotation(foo = anno_horizon(TE_df_cotrans))
library(ComplexHeatmap)
plot_df <- t(TE_df_cotrans)
pdf("cotrans_workflow_HK2.pdf",height=1,width=5.5)
myColor=c("#FDAE61" ,"#FEE090", "#E0F3F8" ,"#74ADD1" ,"#4575B4")
#myColor=c("#4575B4" ,"#74ADD1","#E0F3F8", "#FEE090","#FEE090" ,"#FDAE61" ,"#FDAE61")
Heatmap(plot_df,show_row_dend = FALSE,show_column_names = FALSE,
        column_title = "Bio-sources",column_title_side = "bottom",
        show_column_dend = FALSE,row_names_gp = gpar(fontsize = 8),
        col=myColor,
        heatmap_legend_param = list(
          at = c(-4,-2,0,2,4),
          labels = c(-4,-2,0,2,4),labels_gp = gpar(fontsize = 8),
          title = "TE level", title_position = "topcenter"
        ))
dev.off()

load("data/human_TE_rho.rda")
cor_plot <- human_TE_rho[gene,gene]
col_fun_new_human=colorRamp2(c(0,0.2,0.4,0.6,0.8,1),
                             c("#E1D6D3","#C2B3CB","#9780A1","#9780A1","#9780A1","#60405A"))
tmp_human <- Heatmap(cor_plot,clustering_method_rows = "single",clustering_method_columns="single")
od =  row_order(tmp_human)
human_out = cor_plot[od, od]
human=Heatmap(human_out,row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              col=col_fun_new_human,
              width = ncol(human_out)*unit(3, "mm"), 
              height = nrow(human_out)*unit(3, "mm"),
              heatmap_legend_param = list(
                at = c(0,0.2,0.4,0.6,0.8,1),
                labels = c(0,0.2,0.4,0.6,0.8,1),labels_gp = gpar(fontsize = 8),
                title = "Correlation", title_position = "topcenter"
              ),show_column_names=TRUE,row_names_side = "left",
              cluster_rows = FALSE, cluster_columns = FALSE)
pdf("heatmap_cor.pdf",
    width = 3, height = 2)
draw(human)
dev.off()

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat,diag = TRUE)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

cotran_net <- flattenCorrMatrix(cor_plot)
write.csv(cotran_net,"cotran_net.csv")
