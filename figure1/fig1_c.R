library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(ggplot2)
label_size = 10
legedn_size = 10

setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure1/")

human_all <- read.csv("data/human_QC_summary.csv",row.names=1)
human_anno <- read.csv("data/human_qc_85_filtered_gsm_anno.csv",row.names = 1)
human_anno <- data.frame(t(human_anno))
human_anno_tmp <- merge(human_anno,human_all,by=0)
human_anno_fin <- human_anno_tmp %>% select(Row.names,status)
human_anno_fin$status=sub("keep","Pass",human_anno_fin$status)
human_anno_fin$status=sub("QC filtered out","Fail",human_anno_fin$status)
rownames(human_anno_fin)=human_anno_fin$Row.names

cancer_data <- read.csv("data/cancer_dist.csv",row.names=1)
cancer_data <- as.matrix(cancer_data)
cancer_data <- scale(cancer_data)
cancer_anno_QC <- subset(human_anno_fin, rownames(human_anno_fin) %in% colnames(cancer_data))
cancer_anno_QC <- cancer_anno_QC[order(match(rownames(cancer_anno_QC),colnames(cancer_data))),]

sample_infor =read.csv("data/all_infor_RO1.csv")
cancer_anno_enzyme <- subset(sample_infor, sample_infor$sample %in% colnames(cancer_data))
rownames(cancer_anno_enzyme) = cancer_anno_enzyme$sample
cancer_anno_enzyme <- cancer_anno_enzyme[order(match(cancer_anno_enzyme[,'sample'],colnames(cancer_data))),]

col_fun = colorRamp2(c(-1, 4), c("white", "#984ea3"))
#col_fun = colorRamp2(c(0, 2954124), c("white", "#984ea3"))
human_column_QC = HeatmapAnnotation(QC=cancer_anno_QC$status,annotation_name_gp= gpar(fontsize = 10,fontface = "bold"),
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                                   labels_gp = gpar(fontsize = 8)),
                                    col = list(QC = c("Pass" = "#DFC286", "Fail" = "darkgrey"))) #7DC8CA


human_column_enzyme = HeatmapAnnotation(Digestion=cancer_anno_enzyme$digestion,annotation_name_gp= gpar(fontsize = 10,fontface = "bold"),
                                        annotation_legend_param = list(title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                                       labels_gp = gpar(fontsize = 8)),
                                        col = list(Digestion = c("RNAse I" = "#C07A92", "Mnase" = "#70A3C4","Unknown"="lightgrey")))

human_column_anno = HeatmapAnnotation(QC=cancer_anno_QC$status,Digestion=cancer_anno_enzyme$digestion,
                                      annotation_name_gp= gpar(fontsize = 10,fontface = "bold"),
                                      annotation_legend_param = list(title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                                     labels_gp = gpar(fontsize = 8)),
                                      col = list(QC = c("Pass" = "#DFC286", "Fail" = "darkgrey"),
                                                 Digestion = c("RNAse I" = "#C07A92", "Mnase" = "#70A3C4","Unknown"="lightgrey")))#CFD99D

pdf("heatmap_human_read_length_dist_enzyme_heatscale.pdf",
    width = 5, height = 2.5)

#tiff("heatmap_human_read_length_dist_enzyme.tiff",
    #width = 15, height = 4, res=400,units="in")

draw(Heatmap(cancer_data, cluster_rows = FALSE, use_raster = FALSE, col = col_fun, show_column_names = FALSE, 
        top_annotation = human_column_anno, show_column_dend = FALSE,
        row_title = "Read Length",row_names_gp = gpar(fontsize = 8),row_title_gp = gpar(fontsize = 10)))
dev.off()
