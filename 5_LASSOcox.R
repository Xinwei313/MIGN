

###### Install packages in R 
source("http://bioc.ism.ac.jp/biocLite.R")
biocLite("AnnotationDbi")
install.packages("AnnotationDbi")
install.packages("bit")
biocLite("org.Hs.eg.db")
library(bit)
library(AnnotationDbi)
library(org.Hs.eg.db)
###### Survival package
install.packages("lattice")
library(lattice)


######## plot KM curves
install.packages("reshape2")
library(reshape2)
install.packages("data.table")
library(data.table)
install.packages("zoo")
library("zoo")
install.packages("survminer")
library("survminer")
install.packages("survival")
library("survival")

#### ROC package
install.packages("gplots")
library(gplots)
install.packages("ROCR")
library(ROCR)
install.packages("dplyr")
library("dplyr")
install.packages("lubridate")
library(lubridate)
install.packages("robustbase")
install.packages("caret")
library(caret)
install.packages("pROC")
library(pROC)

install.packages("timeROC")
install.packages("prodlim")
library(prodlim)
library(quantreg)
install.packages("quantreg")
install.packages("polspline")
library("polspline")
library(timeROC)
install.packages("foreign")
library(foreign)

# install.packages("rms")
library(rms)


detach("package:org.Hs.eg.db")
detach("package:foreign")

install.packages("glmnet")
library(glmnet)  # manually load this package by downloading and installing from Toll. 


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


##### (1) TCGA - LGG data 
Survival_LGG=read.csv("TCGA_LGG_Survival data.csv")
Survival_LGG=Survival_LGG[,c(1,2,4,20,75)]
Survival_LGG=na.omit(Survival_LGG)   
rownames(Survival_LGG)=Survival_LGG[,1]
colnames(Survival_LGG)=c('SampleID','OS_event','OS_time','Age','Grade')
Survival_LGG$Grade=substr(Survival_LGG$Grade,2,2)

Gene_LGG=read.csv("TCGA_LGG_GeneExpression.csv")
rownames(Gene_LGG)=Gene_LGG[,1]
colnames(Gene_LGG)=gsub(".", "-", colnames(Gene_LGG), fixed = TRUE)
Gene_LGG=Gene_LGG[,-1]
Gene_LGG=as.matrix(Gene_LGG)

Gene_LGG=Gene_LGG[,intersect(colnames(Gene_LGG),rownames(Survival_LGG))]
Survival_LGG=Survival_LGG[intersect(colnames(Gene_LGG),rownames(Survival_LGG)),]

Gene_LGG=t(Gene_LGG)


###
Q1=apply(Gene_LGG,2,quantile)[2,]
Q2=apply(Gene_LGG,2,quantile)[3,]
Q3=apply(Gene_LGG,2,quantile)[4,]

IQR=Q3-Q1

Gene_LGG=apply(Gene_LGG,1,function(geneExp){
  (geneExp-Q2)/IQR
})

Gene_LGG=na.omit(Gene_LGG)
#Gene_LGG=t(Gene_LGG)  
dim(Gene_LGG)       # 18484   525


Survival_LGG$Grade=as.numeric(Survival_LGG[,5])
Survival=na.omit(Survival_LGG)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


Marker_genes = readRDS("Ligand.rds")

Marker_genes = readRDS("Receptor.rds")

Marker_genes = readRDS("hub_gene.rds")


Survival_Gene=Gene_LGG[t(Marker_genes),rownames(Survival)]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])


Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 525


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


### COX Model
surv=Surv(as.double(time),as.double(status))

set.seed(123456)
cvfit<-cv.glmnet(t(Gene_marker),surv,family="cox",alpha=1,nfolds=10)
dev.new()
plot(cvfit)
A=coef(cvfit, s = "lambda.min")
#A=coef(cvfit, s="lambda.1se")
A=as.numeric(A)
A[(A!=0)]
which=which(A!=0)
sum(A!=0)
Marker_genes[which(A!=0)]


### Ligand
Marker_genes=c("CD24", "EFNB2",  "IGF1", "TNF",  "C3", "HBEGF",  "CXCL11")  # 7 

### Receptor
Marker_genes=c("CHRNA4", "TNFRSF21", "GFRA1", "F11R", "ITGA7", "FXYD6",  "LRRC4",  "P2RY12")  # 8

### HUb genes
Marker_genes3=c("SOX6", "ACVR2B", "DUSP26" , "H2AFY2", "ATCAY",  "CSMD3",  "PAK7", "CUX2", "SHISA7", "MEX3B",  "FAM57B", "MSR1", "ASPM") #13



