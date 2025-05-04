setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure6/")
load("data/cor_TE_exocyst.rda")
load("data/cor_RNA_exocyst.rda")
cor_RNA_exocyst[cor_RNA_exocyst<0]=0
cor_TE_exocyst[cor_TE_exocyst<0]=0
library(circlize)
library(ComplexHeatmap)
library("RColorBrewer")
col_fun_new_human=colorRamp2(c(0,0.2,0.4,0.6,0.8,1),
            c("#45837FFF","#699896FF","#E394BBFF","#DE6DA7FF","#D4419EFF","#D4419EFF"))

#col_fun_new_human=paletteer::paletteer_d("beyonce::X51")
tmp_human <- Heatmap(cor_TE_exocyst,clustering_method_rows = "single",clustering_method_columns="single")
od =  row_order(tmp_human)
human_out = cor_TE_exocyst[od, od]
ha = HeatmapAnnotation(Median_cor = 
                         anno_points(
                           apply(human_out, 1, function(x) { (sum(x) - 1) / 7 }),
                           height = unit(0.5, "cm"),
                           ylim=c(-0,0.6)),
                       gp=gpar(fontsize = 88))

human =Heatmap(human_out,row_names_gp = gpar(fontsize = 6),rect_gp = gpar(type = "none"),
                 col=col_fun_new_human,
               show_heatmap_legend = FALSE,
                width = ncol(human_out)*unit(3, "mm"), 
                height = nrow(human_out)*unit(3, "mm"),
                 heatmap_legend_param = list(
                   at = c(0,0.2,0.4,0.6,0.8,1),
                   labels = c(0,0.2,0.4,0.6,0.8,1),labels_gp = gpar(fontsize = 8),
                   title = "Spearman  ", title_position = "topcenter"
                 ),show_column_names=FALSE,
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(i <= j) {
                     grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                   }
                 })

human_new =Heatmap(human_out,row_names_gp = gpar(fontsize = 6),
                   column_names_gp = gpar(fontsize = 6),
               col=col_fun_new_human,
               show_heatmap_legend = FALSE,
               width = ncol(human_out)*unit(3, "mm"), 
               height = nrow(human_out)*unit(3, "mm"),
               heatmap_legend_param = list(
                 at = c(0,0.2,0.4,0.6,0.8,1),
                 labels = c(0,0.2,0.4,0.6,0.8,1),labels_gp = gpar(fontsize = 8),
                 title = "Spearman  ", title_position = "topcenter"
               ),show_column_names=TRUE,
               cluster_rows = FALSE, cluster_columns = FALSE)
pdf("heatmap_TE_new.pdf",
    width = 5, height = 3)
draw(human_new)
dev.off()

mouse_out = cor_RNA_exocyst[od, od]
ha_m = HeatmapAnnotation(Median_cor = 
                         anno_points(
                           apply(mouse_out, 1, function(x) { (sum(x) - 1) / 7 }),
                           height = unit(0.5, "cm"),
                           ylim=c(0,0.6)),
                       gp=gpar(fontsize = 88))

mouse=Heatmap(mouse_out,column_names_gp = gpar(fontsize = 6),rect_gp = gpar(type = "none"),
        col=col_fun_new_human,
        width = ncol(mouse_out)*unit(3, "mm"), 
        height = nrow(mouse_out)*unit(3, "mm"),
        heatmap_legend_param = list(
          at = c(0,0.2,0.4,0.6,0.8,1),
          labels = c(0,0.2,0.4,0.6,0.8,1),labels_gp = gpar(fontsize = 8),
          title = "Spearman  ", title_position = "topcenter"
        ),show_row_names=FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        })
mouse_new=Heatmap(mouse_out,column_names_gp = gpar(fontsize = 6),
                  row_names_gp = gpar(fontsize = 6),
              col=col_fun_new_human,show_heatmap_legend = FALSE,
              width = ncol(mouse_out)*unit(3, "mm"), 
              height = nrow(mouse_out)*unit(3, "mm"),
              heatmap_legend_param = list(
                at = c(0,0.2,0.4,0.6,0.8,1),
                labels = c(0,0.2,0.4,0.6,0.8,1),labels_gp = gpar(fontsize = 8),
                title = "Spearman  ", title_position = "topcenter"
              ),show_row_names=TRUE,
              cluster_rows = FALSE, cluster_columns = FALSE)

pdf("heatmap_coexpression_new.pdf",
    width = 5, height = 3)
draw(mouse_new)
dev.off()

pdf("heatmap_cor_new.pdf",
width = 5, height = 3)
draw(human + mouse, ht_gap = unit(-55, "mm"))
dev.off()

yarrr::piratepal("all")
