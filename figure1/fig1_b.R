setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure1/")
library(tidyverse)
library(ggrepel)
infor <- read.csv("data/ribobase_7_17_23_manuscript_version.csv")
infor=data.frame(infor)
infor %>% group_by(organism) %>% summarize(dplyr::n())

infor %>% group_by(organism) %>% distinct(study_id) %>% summarize(n())

infor %>% group_by(organism) %>% distinct(cell_line) %>% summarize(n())

###human
tmp1 <- infor %>% filter(organism=="Homo sapiens") %>% group_by(cell_line) %>%
  summarize(biosource_count=n()) %>% arrange(desc(biosource_count))
tmp2 <- infor %>% filter(organism=="Homo sapiens") %>% group_by(cell_line) %>%
  distinct(study_id) %>%
  summarize(cl_study=n()) %>% arrange(desc(cl_study))
human_plot <- merge(tmp1,tmp2, by.x="cell_line",by.y="cell_line")
human_plot_fin <- human_plot %>% mutate(type="Human") %>% 
  arrange(desc(biosource_count)) %>% head(5)
###mouse
tmp1 <- infor %>% filter(organism=="Mus musculus") %>% group_by(cell_line) %>%
  summarize(biosource_count=n()) %>% arrange(desc(biosource_count))
tmp2 <- infor %>% filter(organism=="Mus musculus") %>% group_by(cell_line) %>%
  distinct(study_id) %>%
  summarize(cl_study=n()) %>% arrange(desc(cl_study))
mouse_plot <- merge(tmp1,tmp2, by.x="cell_line",by.y="cell_line")
mouse_plot_fin <- mouse_plot %>% mutate(type="Mouse") %>% 
  arrange(desc(biosource_count)) %>% head(5)
fin_plot <- rbind(human_plot_fin,mouse_plot_fin)
label_size = 10
legedn_size = 10

ggplot(fin_plot,aes(x=cl_study,y=biosource_count,color=type)) + 
  geom_point()+ theme_bw() + ylim(0,300)+ xlim(0,40)+
  geom_text_repel(aes(label =cell_line),size=3) + 
  scale_color_manual(values=c("#984ea3","#4daf4a"))+
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="none",
        legend.text = element_text(size=label_size))+
  labs(y= "Sample numbers",x="Study numbers")
ggsave("figure1_b.pdf", width = 2.4, height = 2.2, units = "in")

