

#=====================================================================================
#
#  Code chunk 1
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
#  Code chunk 2
#
#=====================================================================================


Marker_genes=c("EFNB2","C3","IGF1","TNF", "ASPM", "SOX6", "ACVR2B", "LRRC4", "ATCAY", "FXYD6", 
               "GFRA1", "MSR1", "PAK7", "DUSP26", "TNFRSF21")


Survival_Gene=Gene_CGGALGG[t(Marker_genes),]
Survival=as.matrix(Survival_CGGALGG)

status=as.numeric(Survival[,2])
time=as.numeric(Survival[,3])
Age=as.numeric(Survival[,4])
Grade=as.numeric(Survival[,5])

Genes=rownames(na.omit(Survival_Gene))
Gene_marker=matrix(as.double(as.matrix(Survival_Gene[Genes,])),dim(Survival_Gene)) #58 172


##### univariate Cox #####
p=matrix(0,1,length(Marker_genes))

for (i in 1:length(Marker_genes))
{
  Gene_marker0=Gene_marker[i,]
  
  fit <- coxph(Surv(time, status) ~ (Gene_marker0), method="breslow")
  coeffs <- coef(summary(fit))
  p[i]=coeffs[,5]
}
p


#=====================================================================================
#
#  Code chunk 3
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


#=====================================================================================
#
#  Code chunk 4
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

saveRDS(groups, file = "groupsCGGA.rds")


##############


fit<- survfit(Surv(time/365, status) ~ groups, data = as.data.frame(Survival_Gene))
# Drawing survival curves
dev.new()
ggsurvplot(fit, legend = c(0.2, 0.2), xlim=c(0,7), break.time.by=2)

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
#  Code chunk 5
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
interval=0.872142+c(-1,1)*2*auc.se

#0.8471799 0.8971041


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
interval=0.8986628+c(-1,1)*2*auc.se

#0.8766538 0.9206718


