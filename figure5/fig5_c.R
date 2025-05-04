library(ComplexHeatmap)
library(tidyverse)
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure5/")
glycolysis_list=c("ENO1","ENO2","ENO3","FOXK1",
                  "FOXK2","HK1","HK2","PFKL",
                  "PFKM","PFKP","PGK1","PKM","TPI1","LRRC28")

TE=read.csv("../figure2/data/human_TE_cellline_all_plain.csv",row.names=1)
RNA=read.csv("../figure2/data/human_RNA_cellline_all_plain.csv",row.names=1)

TE_glycolysis=TE[,glycolysis_list]
RNA_glycolysis=RNA[,glycolysis_list]



library(circlize)
tmp=t(TE_glycolysis)
base_mean = rowMeans(tmp)
mat_scaled = t(scale(t(tmp)))
order=data.frame(scale(t(tmp))) %>% arrange(-LRRC28)
mat_scaled=mat_scaled[,rownames(order)]
#Heatmap(mat_scaled)

tmp_RNA=t(RNA_glycolysis)
base_mean_rna = rowMeans(tmp_RNA)
mat_scaled_rna = t(scale(t(tmp_RNA)))
mat_scaled_rna=mat_scaled_rna[,rownames(order)]
Heatmap(mat_scaled_rna)


ha_te = rowAnnotation(foo = anno_block(gp = gpar(fill = "#984ea3"),
                                       labels = "TE", 
                                       labels_gp = gpar(col = "black", fontsize = 10),
                                       ),width = unit(0.5, "cm"))

ha_rna = rowAnnotation(foo = anno_block(gp = gpar(fill = "grey"),
                                        labels = "RNA",
                                        height = unit(10, "mm"),
                                        labels_gp = gpar(col = "black", fontsize = 10)),
                       width = unit(0.5, "cm"))

cluster_columns = FALSE
col_fun = paletteer::paletteer_c("oompaBase::redgreen",n=10)
ht=Heatmap(mat_scaled,row_names_gp = gpar(fontsize = 8),
           cluster_columns = FALSE,
           heatmap_legend_param = list(
             at = c(-4,-2,0,2,4),
             labels = c(-4,-2,0,2,4),labels_gp = gpar(fontsize = 6)),
        cluster_rows = FALSE,right_annotation = ha_te) %v% 
  Heatmap(mat_scaled_rna,row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          cluster_columns = FALSE,
          heatmap_legend_param = list(
            at = c(-4,-2,0,2,4),
            labels = c(-4,-2,0,2,4),labels_gp = gpar(fontsize = 6)),
          cluster_rows = FALSE,show_heatmap_legend = FALSE,
          right_annotation = ha_rna)
pdf("glycolysis.pdf",height=5.5,width=8)
pdf("glycolysis_order.pdf",height=4.5,width=8)
draw(ht)
dev.off()



