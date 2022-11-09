


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


Marker_genes=c("CD24", "EFNB2",  "IGF1", "TNF",  "C3", "HBEGF",  "CXCL11")  # 7 

Marker_genes2=c("CHRNA4", "TNFRSF21", "GFRA1", "F11R", "ITGA7", "FXYD6",  "LRRC4",  "P2RY12")  # 8

Marker_genes3=c("SOX6", "ACVR2B", "DUSP26" , "H2AFY2", "ATCAY",  "CSMD3",  "PAK7", "CUX2", "SHISA7", "MEX3B",  "FAM57B", "MSR1", "ASPM") #13

Marker_genes=unique(c(Marker_genes2, Marker_genes3, "EFNB2","C3","IGF1","TNF"))

#####
#####

Marker_genes=c("EFNB2","C3","IGF1","TNF", "ASPM", "SOX6", "ACVR2B", "LRRC4", "ATCAY", "FXYD6", 
               "GFRA1", "MSR1", "PAK7", "DUSP26", "TNFRSF21")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


Survival_Gene=Gene_LGG[t(Marker_genes),rownames(Survival)]

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])


Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 525


##### #####
dim(Gene_marker)

my.data=t(Gene_marker)
colnames(my.data)=Marker_genes
my.data=as.data.frame(my.data)

my.data$Time=time
my.data$Status=status
my.data$age=Age
my.data$grade=Grade


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


install.packages("My.stepwise")

library(My.stepwise)

?My.stepwise.coxph


my.variable.list = c(Marker_genes,"age","grade")


stepwiseCox = My.stepwise.coxph(Time = "Time", Status = "Status", variable.list = my.variable.list,
                                data = my.data, sle = 0.15, sls = 0.15)


stepwiseCox = My.stepwise.coxph(Time = "Time", Status = "Status", variable.list = my.variable.list, 
                                in.variable = c("EFNB2","C3","IGF1","TNF","age","grade"),
                                data = my.data, sle = 0.15, sls = 0.15)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


##### univariate Cox #####
p=matrix(0,1,length(Marker_genes))

for (i in 1:length(Marker_genes))
{
  Gene_marker0=Gene_marker[i,]
  
  fit <- coxph(Surv(time, status) ~ (Gene_marker0), method="breslow")
  coeffs <- coef(summary(fit))
  p[i]=coeffs[,5]
}
sum(p<0.1)  #58

Marker_genes=Marker_genes[p<0.1]

Gene_marker=Gene_marker[p<0.1,]


##### univariate Cox adjusted with grade #####
p=matrix(0,1,length(Marker_genes))

for (i in 1:length(Marker_genes))
{
  Gene_marker0=Gene_marker[i,]
  Gene_marker0=rbind(Gene_marker0,Grade,Age)
  
  fit <- coxph(Surv(time, status) ~ t(Gene_marker0), method="breslow")
  coeffs <- coef(summary(fit))
  p[i]=coeffs[1,5]
}
sum(p<0.1)  # 149

Marker_genes=Marker_genes[p<0.1]

Gene_marker=Gene_marker[p<0.1,]


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


Gene_marker=rbind(Gene_marker,Age,Grade)
Marker_genes=c(Marker_genes,"Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes


Gene_marker=rbind(Age,Grade)
Marker_genes=c("Age","Grade")
dim(Gene_marker)
rownames(Gene_marker)=Marker_genes


fit <- coxph(Surv(time, status) ~ t(Gene_marker), method="breslow")
coeffs <- coef(summary(fit))

summary(fit)

A=coeffs[,1]

saveRDS(A, file = "trainCoeff_A.rds")


#=====================================================================================
#
#  Code chunk 7 (Forest Plot)
#
#=====================================================================================


multiCox <- summary(fit)

out_multi <- data.frame()
out_multi <- cbind(
  coef = multiCox$coefficients[,"coef"],
  HR = multiCox$conf.int[,"exp(coef)"],
  HR.95L = multiCox$conf.int[,"lower .95"],
  HR.95H = multiCox$conf.int[,"upper .95"],
  pvalue = multiCox$coefficients[,"Pr(>|z|)"]
)


Cell = c("TAM", "TAM", "TAM", "TAM", "TC", "TC", "TC", "TC", "TC", "TC", "TAM", "TC", "TAM", "TC", "TAM", "", "")

Type = c("ligand", "ligand", "ligand", "ligand",  "receptor", "receptor",  "receptor", "receptor",
         "hub gene", "hube gene", "hub gene", "hube gene", "hub gene", "hube gene", "hub gene", "clinical factor", "clinical factor")


out_multi <- as.data.frame(cbind(id=Marker_genes, Cell, Type, out_multi))


out_multi[,4:ncol(out_multi)] <- as.numeric(unlist(out_multi[,4:ncol(out_multi)]))
hz <- paste(round(out_multi$HR,3),
            "(",round(out_multi$HR.95L,3),
            "-",round(out_multi$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Gene",out_multi$id),
                   c(NA,"Cell",out_multi$Cell),
                   c(NA,"Type",out_multi$Type),
                   c(NA,"Coefficient",round(out_multi$coef,3)),
                   c(NA,"P value",round(out_multi$pvalue,3)),
                   c(NA,"Hazard Ratio(95% CI)",hz))


library(forestplot)

dev.new()
forestplot(labeltext=tabletext, 
           graph.pos=3,  #为Pvalue箱线图所在的位置
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,out_multi$HR),
           lower=c(NA,NA,out_multi$HR.95L), #95%置信区间下限
           upper=c(NA,NA,out_multi$HR.95H), #95%置信区间上限
           boxsize=0.3,lwd.ci=2,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=1,      #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0.5, 1,1.5), #横坐标刻度
           lwd.xaxis=1,            #X轴线宽
           lineheight = unit(0.98,"cm"), #固定行高
           graphwidth = unit(.3,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                           "20" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1.5),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


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

saveRDS(groups, file = "groupsTCGA.rds")


##############


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(Survival_Gene))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2), xlim=c(0,5), break.time.by=2)

# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
           xlab='Time (Years)', xlim=c(0,6),break.time.by=2,
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Risk",
           legend.labs = c("Low","High"))


##############


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(Survival_Gene))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2))
# Add risk table
# and change risk table y text colors by strata
ggsurvplot(fit, pval = TRUE, conf.int = FALSE,
           xlab='Time (Years)',
           risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "none",
           legend.title = "Risk",
           legend.labs = c("Low","High"))


########### time dependent ROC 
ROC.DSST<-timeROC(T=time,delta=status,
                  marker=as.numeric(predicted),cause=1,
                  weighting="cox",
                  times=c(1,3,5)*365,ROC=TRUE)
ROC.DSST
dev.new()
plot( ROC.DSST$FP[,1], ROC.DSST$TP[,1],type = "l", lty = 1, pch = NA, col = "blue",ylab="", xlab="");
lines( ROC.DSST$FP[,2], ROC.DSST$TP[,2],type = "l", lty = 1, pch = NA, col = "red"); 
lines( ROC.DSST$FP[,3], ROC.DSST$TP[,3],type = "l", lty = 1, pch = NA, col = "green"); 
lines( c(0,1), c(0,1),type = "l", lty = 1, pch = NA, col = "black"); 
text.legend=c(paste("AUC at 1 years:",as.character(round(ROC.DSST$AUC[1],3))),
              paste("AUC at 3 years:",as.character(round(ROC.DSST$AUC[2],3))),
              paste("AUC at 5 years:",as.character(round(ROC.DSST$AUC[3],3))) )
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)
title(main="Time depedent ROC", ylab="True positive rate", xlab="False positive rate") 
legend("bottomright",pch=c(NA,NA,NA),lty=c(1,1,1),legend=text.legend,col=c("blue","red","green"),bty="n",ncol=1)


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# 90% of the sample, without replacement
n=dim(Gene_marker)[2]
B=1000
aucs=matrix(0,1,B)
for (j in 1:B)
{
  index=sample(c(1:n),size=0.9*n, replace=F)
  Sdata=Gene_marker[,index]
  sta=status[index]
  ti=time[index]
  
  
  predicted=matrix(0,1,dim(Sdata)[2])
  for (i in 1:dim(Sdata)[2])
  {
    predicted[i]=sum(A*Sdata[,i])
  }
  predicted=as.numeric(predicted)
  
  
  
  ROC.DSST<-timeROC(T=ti,delta=sta,
                    marker=as.numeric(predicted),cause=1,
                    weighting="cox",
                    times=c(1,3)*365,ROC=TRUE)
  aucs[j]=ROC.DSST$AUC[1]
}

auc.se=sd(aucs)
interval=0.8800026+c(-1,1)*2*auc.se

#0.8517681 0.9082371


###############
# 90% of the sample, without replacement
n=dim(Gene_marker)[2]
B=1000
aucs=matrix(0,1,B)
for (j in 1:B)
{
  index=sample(c(1:n),size=0.9*n, replace=F)
  Sdata=Gene_marker[,index]
  sta=status[index]
  ti=time[index]
  
  
  predicted=matrix(0,1,dim(Sdata)[2])
  for (i in 1:dim(Sdata)[2])
  {
    predicted[i]=sum(A*Sdata[,i])
  }
  predicted=as.numeric(predicted)
  
  
  
  ROC.DSST<-timeROC(T=ti,delta=sta,
                    marker=as.numeric(predicted),cause=1,
                    weighting="cox",
                    times=c(1,3)*365,ROC=TRUE)
  aucs[j]=ROC.DSST$AUC[2]
}

auc.se=sd(aucs)
interval=0.9029739+c(-1,1)*2*auc.se

#0.8869012 0.9190466




