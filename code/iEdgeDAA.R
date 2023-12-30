#disease gene
library(EnrichmentBrowser)
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)
library(SPIA)



#function

process <- function(x){
  y <- gsub('hsa:','',x)
  y
}
#----------------------------------------
#--------1.GSE19188 differnential analysis--
#lung cancer


data("GSE18842")
exprs_all <- exprs(GSE18842)
# Add the gene symbol
all.eset <- probe2gene(GSE18842)


#lung cancer
data("GSE19188")
exprs_all <- exprs(GSE19188)
# Add the gene symbol
all.eset <- probe2gene(GSE19188)


data("GSE23878")
exprs_all <- exprs(GSE23878)
# Add the gene symbol
all.eset <- probe2gene(GSE23878)


before.norm <- assay(all.eset)
# Gene normalization
all.eset <- normalize(all.eset, norm.method="quantile")
after.norm <- assay(all.eset)

exprs_all1 <- data.frame(after.norm)
table(colData(all.eset)$Group)
colData(all.eset)$GROUP <- ifelse(colData(all.eset)$Group == "d", 1, 0)
normal <- length(which(colData(all.eset)$GROUP == '0'))
tumor <- length(which(colData(all.eset)$GROUP == '1'))

#
library(DrugDiseaseNet)
library(graph)
gg<-keggGlobalGraph()

kegg_genes <- process(nodes(gg))

#write.csv(kegg_genes,'pathway-gene.csv')

left <- c()
right <- c()

for(i in 1:length(edges(gg))){
  cat('i=',i,'\n')
  if(length(edges(gg)[[i]])>0){
    a <- names(edges(gg)[i])
    b <- edges(gg)[[i]]
    for(j in b){
      left <- c(left, a)
      right <- c(right, j)
    }
  }
  
}

edge <- data.frame(left=left, right = right)

#gse18842
pearson <- c()
for(i in 1:length(edge[,1])){
  cat('i=',i,'\n')
  m<- exprs_all1[process(edge[i,1]),c(45:88)]
  m <- t(m)
  n<- exprs_all1[process(edge[i,2]),c(45:88)]
  n <- t(n)
  pearson <- c(pearson,cor(m,n)[1])
  
}

pearson_normal <- c()
for(i in 1:length(edge[,1])){
  cat('i=',i,'\n')
  m<- exprs_all1[process(edge[i,1]),c(1:44)]
  m <- t(m)
  n<- exprs_all1[process(edge[i,2]),c(1:44)]
  n <- t(n)
  pearson_normal <- c(pearson_normal,cor(m,n)[1])
  
}


#gse19188
pearson <- c()
for(i in 1:length(edge[,1])){
  cat('i=',i,'\n')
  m<- exprs_all1[process(edge[i,1]),c(63:153)]
  m <- t(m)
  n<- exprs_all1[process(edge[i,2]),c(63:153)]
  n <- t(n)
  pearson <- c(pearson,cor(m,n)[1])
  
}

pearson_normal <- c()
for(i in 1:length(edge[,1])){
  cat('i=',i,'\n')
  m<- exprs_all1[process(edge[i,1]),c(1:62)]
  m <- t(m)
  n<- exprs_all1[process(edge[i,2]),c(1:62)]
  n <- t(n)
  pearson_normal <- c(pearson_normal,cor(m,n)[1])
  
}

#gse23878
pearson <- c()
for(i in 1:length(edge[,1])){
  cat('i=',i,'\n')
  m<- exprs_all1[process(edge[i,1]),c(20:38)]
  m <- t(m)
  n<- exprs_all1[process(edge[i,2]),c(20:38)]
  n <- t(n)
  pearson <- c(pearson,cor(m,n)[1])
  
}

pearson_normal <- c()
for(i in 1:length(edge[,1])){
  cat('i=',i,'\n')
  m<- exprs_all1[process(edge[i,1]),c(1:19)]
  m <- t(m)
  n<- exprs_all1[process(edge[i,2]),c(1:19)]
  n <- t(n)
  pearson_normal <- c(pearson_normal,cor(m,n)[1])
  
}



edge$pearson <- pearson
edge$pearson_normal <- pearson_normal

edge$pearson_1 <- edge$pearson-edge$pearson_normal

edge1 <- edge[!is.na(edge$pearson_1),]
a = 0.4

edge1_1 <- edge1[which(edge1$pearson_1>a),]
edge1_2 <- edge1[which(edge1$pearson_1 < -a),]

edge1 <- rbind(edge1_1,edge1_2)



#drug-correlation 
library(signatureSearchData)
path <- system.file("extdata", "cmap_instances_02.txt", package="signatureSearchData")
cmap_inst <- read.delim(path, check.names=FALSE)
comp_list <- sampleList(cmap_inst, myby="CMP_CELL")


mas5df <- readRDS("mas5df.rds")

#mas5df <- data.frame(mas5df)

df <- log2(mas5df)



drugs <- c()
scores <- c()
for (i in seq_along(comp_list)) {
  cat('i=',i,'\n')
  
  sample_set <- unlist(comp_list[[i]])
  repcounts <- sapply(comp_list[[i]], length)
  if(repcounts[1]>2 && repcounts[2]>2){
    
    dfsub <- df[, sample_set]
    
    drug_matrix <- dfsub[,1:repcounts[1]]
    control_m <- dfsub[,(repcounts[1]+1):length(sample_set)]
    
    pearson_drug <- c()
    pearson_control <- c()
    for(j in 1:length(edge1[,1])){
      #cat('j=',j,'\n')
      
      if(process(edge1[j,1]) %in% rownames(drug_matrix)){
        m<- drug_matrix[process(edge1[j,1]),]
        m <- t(m)
        m <- as.numeric(m)
      }else{
        m <- rep(NA, repcounts[1])
      }
      
      
      if(process(edge1[j,2]) %in% rownames(drug_matrix)){
        n<- drug_matrix[process(edge1[j,2]),]
        n <- t(n)
        n <- as.numeric(n)
      }else{
        n <- rep(NA, repcounts[1])
      }
      
      pearson_drug <- c(pearson_drug,cor(m,n)[1])
      
      
      if(process(edge1[j,1]) %in% rownames(control_m)){
        m<- control_m[process(edge1[j,1]),]
        m <- t(m)
        m <- as.numeric(m)
      }else{
        m <- rep(NA, repcounts[2])
      }
      
      
      if(process(edge1[j,2]) %in% rownames(control_m)){
        n<- control_m[process(edge1[j,2]),]
        n <- t(n)
        n <- as.numeric(n)
      }else{
        n <- rep(NA, repcounts[2])
      }
      
      pearson_control <- c(pearson_control,cor(m,n)[1])
      
    }
    
    
    edge1$drug1 <- pearson_drug-pearson_control
    
    edge2_1 <- edge1[which(edge1$drug1>a),]
    edge2_2 <- edge1[which(edge1$drug1 < -a),]
    
    edge2 <- rbind(edge2_1,edge2_2)
    
    drugs <- c(drugs,names(comp_list[i]))
    
    score <- sum(-sign(edge2$pearson_1)*sign(edge2$drug1),na.rm = TRUE)
    scores <- c(scores, score)
    
  }
  
  
  
}


edge_drug <- data.frame(drugs = drugs, scores = scores)
process1 <- function(x){
  y <- strsplit(x,'_')[[1]][1]
  y
}

process2 <- function(x){
  y <- strsplit(x,'_')[[1]][2]
  y
}
edge_drug$drug <- sapply(edge_drug$drugs,process1)
edge_drug$cell_line <- sapply(edge_drug$drugs,process2)

#edge_drug <- read.csv('edge_drug_19188_15_10.csv', header = T, stringsAsFactors = F)


write.csv(edge_drug,'edge_drug_31448_03.csv')#more than 2




edge_MCF <- edge_drug[which(edge_drug$cell_line=='MCF7'),]

#edge_MCF <- edge_drug[which(edge_drug$cell_line=='HL60'),]

#edge_MCF <- edge_drug[which(edge_drug$cell_line=='PC3'),]

edge_MCF1 <- edge_MCF[order(edge_MCF$scores, decreasing = T),]

#edge_MCF2 <- edge_MCF1[!duplicated(edge_MCF1$drug),]

rank_all1 <- edge_MCF1

rank_ours <- c()
for(i in 1:length(rank_all1$drug)){
  rank_ours <- c(rank_ours,tolower(rank_all1$drug[i]))
}


library(pROC)
NUM <- 100

CTD <- read.csv('CTD_BC.csv', header = T, stringsAsFactors = F)

drugname <- c()
for(i in 1:length(CTD$Chemical.Name)){
  drugname <- c(drugname,tolower(CTD$Chemical.Name[i]))
}



#rank_ours <- c()
#for(i in 1:length(edge_MCF1$drug)){
#  rank_ours <- c(rank_ours,tolower(edge_MCF1$drug[i]))
#}

library(pROC)

NUM <- 100
drugs <- rank_ours[1:NUM]

y_true <- c()
for(i in 1:NUM){
  if(drugs[i] %in% drugname){
    x <- 1
  }else{
    x <- 0
  }
  y_true <- c(y_true, x)
}

y_pre <- (rank_all1$score[1:NUM]-min(rank_all1$score))/(max(rank_all1$score)-min(rank_all1$score))

#y_pre <- (rank_all1$cs[1:NUM]-min(rank_all1$cs))/(max(rank_all1$cs)-min(rank_all1$cs))


res.roc <- roc(y_true,y_pre)
plot.roc(res.roc,print.auc = T)


#AUPR
library(PRROC)
library(ROCR)
data("ROCR.simple")
scores <- data.frame(y_pre = y_pre, label = y_true)
pr <- pr.curve(scores.class0=scores[scores$label=="1",]$y_pre,
               scores.class1=scores[scores$label=="0",]$y_pre,
               curve=T)
print(pr)

require(ROCR)
pred.class <- as.integer(y_pre > 0.5) #这里假设prediciton>0.5的预测为1，其余为0.
print(cft <- table(pred.class, y_true)) #先prediction value 后label

tp <- cft[2, 2]
tn <- cft[1, 1]
fp <- cft[2, 1]
fn <- cft[1, 2]
print(accuracy <- (tp + tn)/(tp + tn + fp + fn))
## [1] 0.85
print(sensitivity <- tp/(tp + fn))
## [1] 0.8494624
print(specificity <- tn/(tn + fp))
## [1] 0.8504673

print(recall <- tp/(tp+fn))

print(precision <- tp/(tp+fp))
F1 <- 2*recall*precision/(recall+precision)
print(F1)









#AUPR
AUC1 <- c(0.879,
          0.878,
          0.882,
          0.877,
          0.874,
          0.874,
          0.877,
          0.880,
          0.877,
          0.768
          
)

AUC1 <- c(0.6483,0.6632,0.6277,0.6632,0.6815,0.6726,0.6815,0.6726,0.6668,0.6483)

AUC2 <- c(0.6462,0.6306,0.6377,0.6829,0.6547,0.6155,0.6128,0.6472,0.6476,0.6372)
AUC3 <- c(0.6142,0.6771,0.6114,0.6293,0.6763,0.6317,0.5924,0.6671,0.6987,0.638)
AUC4 <- c(0.6339,0.6149,0.6578,0.5974,0.6717,0.6311,0.5899,0.6138,0.5594,0.5876)
AUC <- c(AUC1, AUC2, AUC3, AUC4)


types <- c(rep('5%',10),rep('10%',10),
           rep('15%',10),rep('20%',10))

data <- data.frame(AUC = AUC, types = types)

data$types <- factor(data$types,
                     levels =c('5%','10%',
                               '15%','20%') )

library(ggplot2)
ggplot(data=data, aes(x=types,y=AUC, fill = types)) +
  stat_boxplot(geom = 'errorbar', width = 0.2)+
  geom_boxplot() +
  geom_hline(yintercept=0.7373, colour="red", linetype="dashed", size = 1)+
  xlab('Removal data')+
  scale_y_continuous(breaks = seq(0,1, by = 0.02))+
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
  theme(axis.text.x =element_text(size=12,angle = 45, hjust = 0.5, vjust = 0.5), axis.text.y=element_text(size=12))



#precision
Pre1 <- c(0.82,0.92,0.78,0.86,0.84,0.82,0.84,0.82,0.84,0.82)

Pre1 <- c(0.82,0.86,0.78,0.86,0.84,0.82,0.84,0.82,0.84,0.82)

Pre2 <- c(0.82,0.8,0.84,0.86,0.84,0.8,0.84,0.8,0.76,0.8)
Pre3 <- c(0.8,0.86,0.8,0.78,0.82,0.78,0.76,0.84,0.84,0.8)
Pre4 <- c(0.8,0.82,0.8,0.76,0.78,0.78,0.82,0.84,0.68,0.76)
Precision <- c(Pre1, Pre2, Pre3, Pre4)


types <- c(rep('5%',10),rep('10%',10),
           rep('15%',10),rep('20%',10))

data <- data.frame(Precision = Precision, types = types)

data$types <- factor(data$types,
                     levels =c('5%','10%',
                               '15%','20%') )

library(ggplot2)
ggplot(data=data, aes(x=types,y=Precision, fill = types)) +
  stat_boxplot(geom = 'errorbar', width = 0.2)+
  geom_boxplot() +
  geom_hline(yintercept=0.84, colour="red", linetype="dashed", size = 1)+
  xlab('Removal data')+
  scale_y_continuous(breaks = seq(0,1, by = 0.02))+
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
  theme(axis.text.x =element_text(size=12,angle = 45, hjust = 0.5, vjust = 0.5), axis.text.y=element_text(size=12))
