

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


##### (2) CGGA - LGG data 
Survival_CGGA=read.csv("patient_survival data_CGGA.csv")
Survival_CGGALGG=Survival_CGGA[Survival_CGGA$Grade!=4,]
Survival_CGGALGG=Survival_CGGALGG[,c(1,9,8,5,7)]
Survival_CGGALGG=na.omit(Survival_CGGALGG)
rownames(Survival_CGGALGG)=Survival_CGGALGG[,1]
colnames(Survival_CGGALGG)=c('SampleID','OS_event','OS_time','Age','Grade')

Gene_CGGA=read.csv("Gene expression_CGGA.csv")
Gene_CGGA=as.matrix(Gene_CGGA)
rownames(Gene_CGGA)=Gene_CGGA[,1]
Gene_CGGA=Gene_CGGA[,-1]

Gene_CGGALGG=Gene_CGGA[,intersect(colnames(Gene_CGGA),rownames(Survival_CGGALGG))]            
Survival_CGGALGG=Survival_CGGALGG[intersect(colnames(Gene_CGGA),rownames(Survival_CGGALGG)),] 

Gene_CGGALGG=t(Gene_CGGALGG)


###
sample = rownames(Gene_CGGALGG)
Gene_CGGALGG <- as.data.frame(Gene_CGGALGG)
Gene_CGGALGG <- sapply(Gene_CGGALGG, as.numeric)
rownames(Gene_CGGALGG) <- sample


###
Q1=apply(Gene_CGGALGG,2,quantile)[2,]
Q2=apply(Gene_CGGALGG,2,quantile)[3,]
Q3=apply(Gene_CGGALGG,2,quantile)[4,]

IQR=Q3-Q1

Gene_CGGALGG=apply(Gene_CGGALGG,1,function(geneExp){
  (geneExp-Q2)/IQR
})

Gene_CGGALGG=na.omit(Gene_CGGALGG)
#Gene_CGGALGG=t(Gene_CGGALGG)
dim(Gene_CGGALGG)       # 19151   172


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


##### my signature
Marker_genes=c("EFNB2","C3","IGF1","TNF", "ASPM", "SOX6", "ACVR2B", "LRRC4", "ATCAY", "FXYD6", 
               "GFRA1", "MSR1", "PAK7", "DUSP26", "TNFRSF21")


Survival_Gene=Gene_LGG[t(Marker_genes),rownames(Survival)]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 525

###
Gene_marker=rbind(Gene_marker,Age,Grade)
Marker_genes=c(Marker_genes,"Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
coeffs <- coef(summary(fit))
A=coeffs[,1]


#####
#####


Survival_Gene=Gene_CGGALGG[t(Marker_genes),]
Survival=as.matrix(Survival_CGGALGG)

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 172


###
Gene_marker=rbind(Gene_marker,Age,Grade)
Marker_genes=c(Marker_genes,"Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}
###
groups=matrix(0,1,length(status))
groups[predicted<=median(predicted)]=1
groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)


troc1=timeROC(T=time,delta=status,
              marker=as.numeric(predicted),cause=1,
              weighting="cox",
              times=c(1,3,5)*365,ROC=TRUE)

saveRDS(troc1,"troc1.rds")


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Xiaoqiang
Marker_genes=c('ANPEP','DPP4','PRRG1','GPNMB','TMEM26','PXDN','CDH6','SCN3A','SEMA6B','CCDC37',
               'FANCA','NETO2')


Survival_Gene=Gene_LGG[t(Marker_genes),rownames(Survival)]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 525

###
Gene_marker=rbind(Gene_marker,Age,Grade)
Marker_genes=c(Marker_genes,"Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
coeffs <- coef(summary(fit))
A=coeffs[,1]


#####
#####


Survival_Gene=Gene_CGGALGG[t(Marker_genes),]
Survival=as.matrix(Survival_CGGALGG)

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 172


###
Gene_marker=rbind(Gene_marker,Age,Grade)
Marker_genes=c(Marker_genes,"Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}
###
groups=matrix(0,1,length(status))
groups[predicted<=median(predicted)]=1
groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)


troc2=timeROC(T=time,delta=status,
              marker=as.numeric(predicted),cause=1,
              weighting="cox",
              times=c(1,3,5)*365,ROC=TRUE)

saveRDS(troc2,"troc2.rds")


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Novel Immune-Related Gene Signature for Risk Stratification and Prognosis of Survival in Lower-Grade Glioma
Marker_genes=c("CANX","HSPA1B","KLRC2","PSMC6","RFXAP","TAP1")


Survival_Gene=Gene_LGG[t(Marker_genes),rownames(Survival)]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 525

###
Gene_marker=rbind(Gene_marker,Age,Grade)
Marker_genes=c(Marker_genes,"Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
coeffs <- coef(summary(fit))
A=coeffs[,1]


#####
#####


Survival_Gene=Gene_CGGALGG[t(Marker_genes),]
Survival=as.matrix(Survival_CGGALGG)

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 172


###
Gene_marker=rbind(Gene_marker,Age,Grade)
Marker_genes=c(Marker_genes,"Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}
###
groups=matrix(0,1,length(status))
groups[predicted<=median(predicted)]=1
groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)


troc3=timeROC(T=time,delta=status,
              marker=as.numeric(predicted),cause=1,
              weighting="cox",
              times=c(1,3,5)*365,ROC=TRUE)

saveRDS(troc3,"troc3.rds")


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


Marker_genes=c("CANX","HSPA1B","KLRC2","PSMC6","RFXAP","TAP1")


Survival_Gene=Gene_LGG[t(Marker_genes),rownames(Survival)]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 525

###
Gene_marker=rbind(Age,Grade)
Marker_genes=c("Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
coeffs <- coef(summary(fit))
A=coeffs[,1]


#####
#####


Survival_Gene=Gene_CGGALGG[t(Marker_genes),]
Survival=as.matrix(Survival_CGGALGG)

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 172


###
Gene_marker=rbind(Age,Grade)
Marker_genes=c("Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes

########  ROC for Training set
predicted=matrix(0,1,dim(Gene_marker)[2])
for (i in 1:dim(Gene_marker)[2])
{
  predicted[i]=sum(A*Gene_marker[,i])
}
###
groups=matrix(0,1,length(status))
groups[predicted<=median(predicted)]=1
groups[predicted>median(predicted)]=2
groups=t(groups)
groups=as.numeric(groups)


troc4=timeROC(T=time,delta=status,
              marker=as.numeric(predicted),cause=1,
              weighting="cox",
              times=c(1,3,5)*365,ROC=TRUE)

saveRDS(troc4,"troc4.rds")


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# 
dev.new()
plot( troc1$FP[,1], troc1$TP[,1],type = "l", lty = 1, pch = NA, col = "red",ylab="", xlab="");
lines( troc2$FP[,1], troc2$TP[,1],type = "l", lty = 1, pch = NA, col = "green"); 
lines( troc3$FP[,1], troc3$TP[,1],type = "l", lty = 1, pch = NA, col = "blue"); 
lines( troc4$FP[,1], troc4$TP[,1],type = "l", lty = 1, pch = NA, col = "black"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
text.legend=c(paste("AUC of signature 1:",as.character(round(troc1$AUC[1],3))),paste("AUC of signature 2:",as.character(round(troc2$AUC[1],3))),
              paste("AUC of signature 3:",as.character(round(troc3$AUC[1],3))),paste("AUC of signature 4:",as.character(round(troc4$AUC[1],3))))
legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("blue","red","black"),bty="n",ncol=1)
title(main="1-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("red","green","blue","black"),bty="n",ncol=1)
# 


# 
dev.new()
plot( troc1$FP[,2], troc1$TP[,2],type = "l", lty = 1, pch = NA, col = "red",ylab="", xlab="");
lines( troc2$FP[,2], troc2$TP[,2],type = "l", lty = 1, pch = NA, col = "green"); 
lines( troc3$FP[,2], troc3$TP[,2],type = "l", lty = 1, pch = NA, col = "blue"); 
lines( troc4$FP[,2], troc4$TP[,2],type = "l", lty = 1, pch = NA, col = "black"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
text.legend=c(paste("AUC of signature 1:",as.character(round(troc1$AUC[2],3))),paste("AUC of signature 2:",as.character(round(troc2$AUC[2],3))),
              paste("AUC of signature 3:",as.character(round(troc3$AUC[2],3))),paste("AUC of signature 4:",as.character(round(troc4$AUC[2],3))))
legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("blue","red","black"),bty="n",ncol=1)
title(main="3-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("red","green","blue","black"),bty="n",ncol=1)
# 


# 
dev.new()
plot( troc1$FP[,3], troc1$TP[,3],type = "l", lty = 1, pch = NA, col = "red",ylab="", xlab="");
lines( troc2$FP[,3], troc2$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
lines( troc3$FP[,3], troc3$TP[,3],type = "l", lty = 1, pch = NA, col = "blue"); 
lines( troc4$FP[,3], troc4$TP[,3],type = "l", lty = 1, pch = NA, col = "black"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "gray"); 
text.legend=c(paste("AUC of signature 1:",as.character(round(troc1$AUC[3],3))),paste("AUC of signature 2:",as.character(round(troc2$AUC[3],3))),
              paste("AUC of signature 3:",as.character(round(troc3$AUC[3],3))),paste("AUC of signature 4:",as.character(round(troc4$AUC[3],3))))
legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("blue","red","black"),bty="n",ncol=1)
title(main="5-year survival prediction ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA,NA),lty=c(1,1,1,1),legend=text.legend,col=c("red","green","blue","black"),bty="n",ncol=1)
# 

