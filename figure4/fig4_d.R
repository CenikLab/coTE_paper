library(tidyverse)
library(ROCit)
library(GO.db)
library(stringr)

label_size = 10
legedn_size = 10

setwd("/Users/yliu5/Downloads/final_submit/figure_code/figure4/")
human_data=read.csv("data/all_TE_better_auroc_mean_human.csv",row.names=1)

GO_all <- as.list(GOTERM)
human_data$anno <-sapply(human_data$term, function(x) {
  GO_all[[x]]@Term
})
#write.csv(human_data,"human_data_anno.csv")
###cor all vs. term numbers

GO="GO:0004709" #MPK kinase activity
my.term <- GO_all$`GO:0004709`@Term
load("../figure3/data/human_TE_rho.rda")
load("data/human_RNA_rho_new.rda")
load("data/annotations_12_4_22_human.rda")

ranking <- function(human_TE_rho,GO,type){
  network=human_TE_rho
  #if(type=="TE"){
  #rownames(annotations_12_4_22_human)=gsub("-",".",rownames(annotations_12_4_22_human))}
  rownames(annotations_12_4_22_human)=gsub("-",".",rownames(annotations_12_4_22_human))
  anno=annotations_12_4_22_human[,GO,drop=FALSE]
  genes.labels=as.matrix(anno)
  network=abs(network)
  ord <- order(rownames(network))
  network <- network[ord, ord]
  ord <- order(rownames(genes.labels))
  genes.labels <- as.matrix(genes.labels[ord, ])
  match.lab <- match(rownames(genes.labels), rownames(network))
  filt.lab <- !is.na(match.lab)
  filt.net <- match.lab[filt.lab]
  network <- network[filt.net, filt.net]
  genes.labels <- as.matrix(genes.labels[filt.lab, ])
  
  l <- dim(genes.labels)[2]
  g <- dim(genes.labels)[1]
  ab <- which(genes.labels != 0, arr.ind = TRUE)
  n <- length(ab[, 1])
  ###mask one gene each time
  j=1
  test.genes.labels <- matrix(nrow=g)
  d <- which(ab[, 2] == j)  # Which indices the genes are in this particular GO group
  t <- length(d)  # Total number of genes in the GO group
  test.genes.labels=matrix(genes.labels[,j], nrow = g, ncol = t)
  ab_tmp=which(genes.labels[,j,drop=FALSE] != 0, arr.ind = TRUE)
  for (i in 1:t) {
    c <- 1 + (i - 1)  # GO group to look at (ie column)
    test.genes.labels[ab_tmp[,'row'], c][i] <- 0
    #tmp[ab_tmp[d], c][i] <- 0
  }
  ###known filter
  #filter <- matrix(genes.labels, nrow = g, ncol = n * l)
  #filter <- matrix(nrow=g)
  d <- which(ab[, 2] == j)  # Which indices the genes are in this particular GO group
  t <- length(d)  # Total number of genes in the GO group
  filter=matrix(genes.labels[,j], nrow = g, ncol = t)
  
  sumin <- crossprod(network, test.genes.labels)
  
  sumall <- matrix(apply(network, 2, sum), ncol = dim(sumin)[2], nrow = dim(sumin)[1])
  
  predicts <- sumin/sumall
  
  nans <- which(test.genes.labels == 1, arr.ind = TRUE)
  predicts[nans] <- NA
  
  predicts <- apply(abs(predicts), 2, rank, na.last = "keep", ties.method = "average")
  
  negatives <- which(filter == 0, arr.ind = TRUE)
  positives <- which(filter == 1, arr.ind = TRUE)
  predicts[negatives] <- 0
  predicts[ab_tmp[,'row'],]
  #predicts_fin=data.frame(x1=sample(x=c(1:g), size=g),row.names = row.names(genes.labels))
  predicts_fin=data.frame(x1=rep(0,g),row.names = row.names(genes.labels))
  for (i in 1:t) {
    c <- 1 + (i - 1)  # GO group to look at (ie column)
    predicts_fin[ab_tmp[i,'row'],] <- predicts[ab_tmp[,'row'], c][i]
    #tmp[ab_tmp[d], c][i] <- 0
  }
  random_tmp=c(1:g)
  random=random_tmp[!random_tmp %in% predicts_fin[ab_tmp[,'row'],]]
  predicts_fin$x1[predicts_fin$x1==0] <- sample(x=random, size=length(random))
  return(predicts_fin)
}
TE_predict=ranking(human_TE_rho,GO,"TE")
RNA_predict=ranking(human_RNA_rho_new,GO,"RNA")

class_ranking <- function(human_TE_rho,type){
    network=human_TE_rho
    rownames(annotations_12_4_22_human)=gsub("-",".",rownames(annotations_12_4_22_human))
    anno=annotations_12_4_22_human[,GO,drop=FALSE]
    genes.labels=as.matrix(anno)
    network=abs(network)
    ord <- order(rownames(network))
    network <- network[ord, ord]
    ord <- order(rownames(genes.labels))
    genes.labels <- as.matrix(genes.labels[ord, ])
    match.lab <- match(rownames(genes.labels), rownames(network))
    filt.lab <- !is.na(match.lab)
    filt.net <- match.lab[filt.lab]
    network <- network[filt.net, filt.net]
    genes.labels <- as.matrix(genes.labels[filt.lab, ])
    return(genes.labels[,1])
}

class_TE <- class_ranking(human_TE_rho,"TE")
class_RNA <- class_ranking(human_RNA_rho_new,"RNA")
human_TE <- rocit(score = TE_predict[,"x1"], class = class_TE) 
human_RNA <- rocit(score = RNA_predict[,"x1"], class = class_RNA) 
roc_te <- summary(human_TE)
roc_rna <- summary(human_RNA)
pdf("HuMAP2_GO:0004709_fin.pdf",
    width =3.04, height = 3.52)
plot(human_TE, legend=F,col="#984ea3",lwd=4,YIndex = FALSE)
lines(human_RNA$TPR~human_RNA$FPR,col="grey",lwd=3)
text(0.7,0.3,paste("TE:",round(as.numeric(sub('Area under curve: ','',roc_te[4,1])),2),sep=""),col="#984ea3",cex=0.8)
text(0.7,0.1,paste("RNA:",round(as.numeric(sub('Area under curve: ','',roc_rna[4,1])),2),sep=""),col="grey",cex=0.8)
title(main="MAP kinase activity")
dev.off()

gene_list=rownames(annotations_12_4_22_human[annotations_12_4_22_human$`GO:0004709`==1,,drop=FALSE])
GO_MPK_TE=human_TE_rho[gene_list,gene_list]
GO_MPK_RNA=human_RNA_rho_new[gene_list,gene_list]
save(GO_MPK_TE,file="GO_MPK_TE.rda")
save(GO_MPK_RNA,file="GO_MPK_RNA.rda")
