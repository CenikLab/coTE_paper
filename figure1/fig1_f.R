setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure1/")
setwd("/Users/yliu5/Library/CloudStorage/Box-Box/Canlab/publication_2023_version/Figure/Figure2/summary_QC/")
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(extrafont)

label_size = 10
legedn_size = 10

all_summary <- read.csv("data/summary_new.csv")
all_plot <- all_summary %>% pivot_longer(c(Human,Mouse),names_to="Type",values_to = "Counts") 

all_plot %>% na.omit() %>%
  group_by(Type, cum_tot) %>% 
  mutate(cum = cumsum(Counts)) %>% ggplot(aes(x=cum_tot, y=cum, fill =Type,label = Counts)) + 
  geom_col(data = . %>% filter( QC=="Pass"), position = position_dodge(width = 0.7), alpha = 1,width=0.6) +
  geom_col(data = . %>% filter( QC=="Fail"), position = position_dodge(width = 0.7), alpha = 0.4,width=0.6) +
  geom_text(data = . %>% filter( QC=="Pass"), position = position_dodge(width = 0.9),hjust = 1,size=3) +
  geom_text(data = . %>% filter( QC=="Fail"), position = position_dodge(width = 0.9),hjust = 1,size=3) +
  geom_tile(aes(y=NA_integer_, alpha = factor(QC))) + 
  scale_alpha_manual(values = c(0.4,1)) +
  theme_bw() +  scale_fill_manual(values=c("#984ea3","#4daf4a")) +
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),legend.position="left",
        legend.text = element_text(size=label_size))+labs(y= "Number of samples",x="")
ggsave("all_data_summary_new.pdf", width = 2.4, height = 2, units = "in")
