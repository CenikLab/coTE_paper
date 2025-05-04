library(tidyverse)
setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure2/")
homo <- read.csv("data/shared_gene_v2.csv",row.names=1)

human_TE <- read.csv("data/human_TE_cellline_all_plain.csv",row.names=1)
mouse_TE <- read.csv("data/mouse_TE_cellline_all_plain.csv",row.names=1)

colnames(human_TE) <- tolower(colnames(human_TE))
colnames(mouse_TE) <- tolower(colnames(mouse_TE))

human_TE_homo <- human_TE[,homo$share_gene]
mouse_TE_homo <- mouse_TE[,homo$share_gene]

label_size = 10
legedn_size = 10

library(ggplot2)
library(umap)
library(dplyr)

combined_data <- rbind(human_TE_homo, mouse_TE_homo)
species <- factor(c(rep("human", nrow(human_TE_homo)), rep("mouse", nrow(mouse_TE_homo))))
umap_result <- umap(combined_data)

# Convert UMAP results to a data frame
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$species <- species

# Plot the UMAP results with colors indicating different species
library(randomcoloR)
umap_df=read.csv("data/umap_df.csv",row.names=1)
umap_df_list <- umap_df %>% group_by(position) %>% summarise(count=n()) %>%
  filter(count>=5) %>% 
  filter(position!="unknown")

umap_df_plot <- umap_df %>% filter(position %in% umap_df_list$position)
newCols <- distinctColorPalette((length(unique(umap_df_plot$position))))
library(ggplot2)
ggplot(umap_df_plot, aes(x = UMAP1, y = UMAP2, color = position, shape = species)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = newCols) +
  theme_bw()+theme(axis.text=element_text(size=label_size),
                   axis.title = element_text(size=label_size),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(), panel.grid.major = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.text = element_text(size=label_size)) +
                   #legend.position = "none")+
  labs(color = "position")
ggsave("human_mouse_uMPA_TE_plot.pdf", width = 2.4, height = 2.2, units = "in")

write.csv(umap_df,"umap_df.csv")

### cancer + specise
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cancer, shape = species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP of Combined Human and Mouse Data", color = "position")

