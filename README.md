# iEdgePathDDA: integrated edge information and pathway topology for drug-disease associations

iEdgeDDA contains the basic R functions and sample data for running the iEdgeDDA algorithm. After installing and loading the package, users will be able to explore the framework of iEdgeDDA


# More about iEdgePathDDA
The iEdgeDDA is a novel tool used to identify the potentail candidate drugs for cancers based on edge information and pathway topology.

# Getting started

## Step 1. Pre run the method for installation

You should ensure that you have the necessary system dependencies configured.

For Windows (8.1 / 10 / 11): Rtools should be installed to the system path.

The latest base R is recommended. The compatibility of the earlier version (v4.0.x) is under evaluation.
We use R version is [64-bit] d:\Program Files\R\R-4.1.2

## Step 2. Install the package
The dependency `EnrichmentBrowser`, `KEGGdzPathwaysGEO`, `KEGGandMetacoreDzPathwaysGEO`, `SPIA` and `signatureSearchData` are unavailable on the CRAN but available on [BioConductor](https://www.bioconductor.org/). So we need to install the BiocManager manually. 

``` r
if (!"BiocManager" %in% as.data.frame(installed.packages())$Package)
  install.packages("BiocManager")
BiocManager::install(c("EnrichmentBrowser", "KEGGdzPathwaysGEO","KEGGandMetacoreDzPathwaysGEO","SPIA","signatureSearchData"))
```
Then you can install the development version of iEdgeDDA from [GitHub](https://github.com/) with:

``` r
if (!"devtools" %in% as.data.frame(installed.packages())$Package)
  install.packages("devtools")
devtools::install_github("eshinesimida/iEdgeDDA")

```
## Examples

Below is a basic example that shows how to obtain candidate drugs of lung cancer:

``` r
#Load require package
library(EnrichmentBrowser)
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)
library(SPIA)



#function

process <- function(x){
  y <- gsub('hsa:','',x)
  y
}
#---------load the data of lung cancer
data("GSE18842")
exprs_all <- exprs(GSE18842)
# Add the gene symbol
all.eset <- probe2gene(GSE18842)
before.norm <- assay(all.eset)
# Gene normalization
all.eset <- normalize(all.eset, norm.method="quantile")
after.norm <- assay(all.eset)

exprs_all1 <- data.frame(after.norm)
table(colData(all.eset)$Group)
colData(all.eset)$GROUP <- ifelse(colData(all.eset)$Group == "d", 1, 0)
normal <- length(which(colData(all.eset)$GROUP == '0'))
tumor <- length(which(colData(all.eset)$GROUP == '1'))


#2.load pathway information from kegg database
library(DrugDiseaseNet)
library(graph)
gg<-keggGlobalGraph()

kegg_genes <- process(nodes(gg))

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

#3. get edge information from pathways
edge <- data.frame(left=left, right = right)

#4. get disease-related edges (the changes in gene interaction)
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

edge$pearson <- pearson
edge$pearson_normal <- pearson_normal
edge$pearson_1 <- edge$pearson-edge$pearson_normal
edge1 <- edge[!is.na(edge$pearson_1),]
a = 0
edge1_1 <- edge1[which(edge1$pearson_1>a),]
edge1_2 <- edge1[which(edge1$pearson_1 < -a),]
edge1 <- rbind(edge1_1,edge1_2)



#5. Calculate the inhibition score between drug-induced edges and disease-related edges
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
edge_MCF <- edge_drug[which(edge_drug$cell_line=='MCF7'),]
edge_MCF1 <- edge_MCF[order(edge_MCF$scores, decreasing = T),]
rank_all1 <- edge_MCF1

