setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure6/")
#summary_df=read.csv("summary_allterms_fin.csv")
summary_df=read.csv("data/summary_allterms_string_5_27_24.csv")
label_size=8
summary_df %>%
  mutate(label = factor(label, 
                        levels=c("s1", "s2","s3","s4","s6", "s7", 
                                 "s8", "s9",
                                 "s5"))) %>% 
  mutate(group = as.factor(group),
         group = fct_relevel(group, "string", "GO_human", "database")) %>% 
  ggplot(aes(x=label,y=percentage,
             fill=group))+geom_col(width=0.5,position = "dodge2")+
  theme_bw() + scale_fill_manual(values=c("#F27983","#8FBDD9","grey")) + 
  theme(axis.text=element_text(size=label_size),
        axis.title = element_text(size=label_size),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        legend.position="top",
        legend.text = element_text(size=label_size))+ labs(y= "Percentage",x="")
ggsave("net_AUROC_string_new_one_level_5_27_24.pdf",width = 3 , height = 2.5, units = "in")

library(rcompanion)
tmp=read.csv("data/chi_sq_test_all.csv",row.names=1)
data=tmp[,c(1,3)]
chi_square_test <- function(observed, all, total_string, total_all) {
  contingency_table <- matrix(c(observed, total_string - observed, all, total_all - all), nrow = 2)
  test <- chisq.test(contingency_table)
  
  # Calculate effect size (Phi coefficient)
  phi <- sqrt(test$statistic / sum(contingency_table))
  
  # Return chi-square statistic, p-value, and effect size
  return(c(test$statistic, test$p.value, phi))
}
results <- data.frame(Type = character(), Chi2 = numeric(), p.value = numeric(), EffectSize = numeric())
total_string <- sum(data$STRING)
total_all <- sum(data$ALL)
for (i in 1:nrow(data)) {
  result <- chi_square_test(data$STRING[i], data$ALL[i], total_string, total_all)
  results <- rbind(results, c(data$Group[i], result))
}

colnames(results) <- c("Chi2", "p.value", "EffectSize")
rownames(results) = rownames(tmp)
write.csv(results,"string_chi_p_effect.csv")

data=tmp[,c(2,3)]
chi_square_test <- function(observed, all, total_string, total_all) {
  contingency_table <- matrix(c(observed, total_string - observed, all, total_all - all), nrow = 2)
  test <- chisq.test(contingency_table)
  
  # Calculate effect size (Phi coefficient)
  phi <- sqrt(test$statistic / sum(contingency_table))
  
  # Return chi-square statistic, p-value, and effect size
  return(c(test$statistic, test$p.value, phi))
}

results <- data.frame(Type = character(), Chi2 = numeric(), p.value = numeric(), EffectSize = numeric())
total_GO <- sum(data$GO)
total_all <- sum(data$ALL)
for (i in 1:nrow(data)) {
  result <- chi_square_test(data$GO[i], data$ALL[i], total_GO, total_all)
  results <- rbind(results, c(data$Group[i], result))
}
colnames(results) <- c("Chi2", "p.value", "EffectSize")
rownames(results) = rownames(tmp)
write.csv(results,"GO_chi_p_effect.csv")

data=tmp[,c(1,2)]

chi_square_test <- function(observed, all, total_string, total_all) {
  contingency_table <- matrix(c(observed, total_string - observed, all, total_all - all), nrow = 2)
  test <- chisq.test(contingency_table)
  
  # Calculate effect size (Phi coefficient)
  phi <- sqrt(test$statistic / sum(contingency_table))
  
  # Return chi-square statistic, p-value, and effect size
  return(c(test$statistic, test$p.value, phi))
}
results <- data.frame(Type = character(), Chi2 = numeric(), p.value = numeric(), EffectSize = numeric())
total_string <- sum(data$STRING)
total_GO <- sum(data$GO)
for (i in 1:nrow(data)) {
  result <- chi_square_test(data$STRING[i], data$GO[i], total_string, total_GO)
  results <- rbind(results, c(data$Group[i], result))
}

colnames(results) <- c("Chi2", "p.value", "EffectSize")
rownames(results) = rownames(tmp)
write.csv(results,"string_go_chi_p_effect.csv")
