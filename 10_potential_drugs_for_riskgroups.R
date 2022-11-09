


install.packages("oncoPredict")
install.packages("car")

BiocManager::install("GenomicFeatures")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("sva")

library(sva)
library(car)

rm(list = ls())  ## Ä§»Ã²Ù×÷£¬Ò»¼üÇå¿Õ~
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./oncoPredict_DataFiles/DataFiles/Training Data/'


GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))

GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 
GDSC2_Res[1:4,1:4]


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

dim(Gene_LGG)
testExpr <- Gene_LGG


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

Gene_CGGALGG=t(Gene_CGGALGG)

dim(Gene_CGGALGG)
testExpr <- Gene_CGGALGG

#####
#####

Survival_CGGA=read.csv("patient_survival data_CGGA.csv")
Survival_CGGA=Survival_CGGA[Survival_CGGA$Grade!=4,]
Survival_CGGA=Survival_CGGA[,c(1,9,8,5,2,7)]
Survival_CGGA=na.omit(Survival_CGGA)
rownames(Survival_CGGA)=Survival_CGGA[,1]
colnames(Survival_CGGA)=c('SampleID','OS_event','OS_time','Age','Subtype','Grade')

#Survival_CGGA=Survival_CGGA[Survival_CGGA$Subtype=="Proneural",]

Gene_CGGA=read.csv("Gene expression_CGGA.csv")
Gene_CGGA=as.matrix(Gene_CGGA)
rownames(Gene_CGGA)=Gene_CGGA[,1]
Gene_CGGA=Gene_CGGA[,-1]

Gene_CGGA=Gene_CGGA[,intersect(colnames(Gene_CGGA),rownames(Survival_CGGA))]            
Survival_CGGA=Survival_CGGA[intersect(colnames(Gene_CGGA),rownames(Survival_CGGA)),] 

subtype=as.factor(Survival_CGGA$Subtype)  # 99
df=t(Gene_CGGA)  
dim(df)     # 99 21496

###
sample = rownames(df)
df <- as.data.frame(df)
df <- sapply(df, as.numeric)
rownames(df) <- sample


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


## read the result
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]


rownames(testPtype)=testPtype[,1]
testPtype=testPtype[,-1]
dim(testPtype)


##### #####
TCGAsample=rownames(Survival)
testPtype=testPtype[TCGAsample,]

groups <- readRDS("groupsTCGA.rds")


##### #####
groups <- readRDS("groupsCGGA.rds")


##### #####
groups = Survival_CGGA$Subtype
groups = as.factor(groups)


##### 
#####

p=matrix(0,1,ncol(testPtype))

for (i in 1:ncol(testPtype))
{
  testPtype0=testPtype[,i]
  
  S1=testPtype0[groups==1]
  S2=testPtype0[groups==2]
  
  Wtest=wilcox.test(S1,S2) 
  
  p[i]=Wtest$p.value
}
sum(p<0.01)  #58

drugCGGA=colnames(testPtype)[p<0.01]  # 46
drugTCGA=colnames(testPtype)[p<0.01]  # 88


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


###TCGA
#igf1r=c("Linsitinib_1510", "NVP-ADW742_1932", "IGF1R_3801_1738", "GSK1904529A_1093"  )

igf1r=c("Linsitinib_1510", "NVP-ADW742_1932", "GSK1904529A_1093"  )


#pi3k=c("GNE-317_1926", "AZD6482_2169", "Dactolisib_1057", "PF-4708671_1129", "Uprosertib_2106", "AT13148_2170", "Pictilisib_1058",
#       "CZC24832_1615", "Buparlisib_1873", "Afuresertib_1912", "Ipatasertib_1924", "LJI308_2107", "MK-2206_1053", "AZD8055_1059" )

pi3k=c("GNE-317_1926", "AZD6482_2169", "Dactolisib_1057", "PF-4708671_1129", "Uprosertib_2106", "Pictilisib_1058",
       "CZC24832_1615", "Buparlisib_1873", "Afuresertib_1912", "Ipatasertib_1924", "MK-2206_1053", "AZD8055_1059" )

intersect(pi3k, drugTCGA)



###CGGA
#igf1r=c("GSK1904529A_1093", "Linsitinib_1510",  "NVP-ADW742_1932", "BMS-754807_2171")

igf1r=c("Linsitinib_1510", "NVP-ADW742_1932", "GSK1904529A_1093"  )


pi3k=c("PF-4708671_1129", "Alpelisib_1560", "Taselisib_1561", "OSI-027_1594", "CZC24832_1615")


intersect(igf1r, drugTCGA)
# "GSK1904529A_1093" "Linsitinib_1510"  "NVP-ADW742_1932" 


intersect(pi3k, drugCGGA)
#"PF-4708671_1129" "CZC24832_1615"

predicted=testPtype[,"Linsitinib_1510"]
predicted=testPtype[,"NVP-ADW742_1932"]
predicted=testPtype[,"GSK1904529A_1093"]

S1=predicted[groups==1]
S2=predicted[groups==2]

S1=predicted[groups=="Mesenchymal"]
S2=predicted[groups=="Classical"]
S3=predicted[groups=="Neural"]
S4=predicted[groups=="Proneural"]


Wtest=wilcox.test(S1,S2)  # alternative = c("two.sided", "less", "greater")  "two.sided" (default),
pvalue=Wtest$p.value
pvalue

dataset1 <- data.frame(value = c(S1,S2), 
                       group = factor(rep(c("Low-risk", "High-risk"), 
                                          times = c(length(S1),length(S2)))))



names(dataset1) <- c('Drug sensitivity score', 'Risk group')


p1 <- ggboxplot(dataset1, x='Risk group', y='Drug sensitivity score', fill='Risk group', palette='npg',
                add='none') 


dev.new()
p1+stat_compare_means(method = "wilcox.test", label.x = 0.56)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


dataset1 <- data.frame(value = c(S1,S2,S3,S4), 
                       group = factor(rep(c("Mesenchymal", "Classical", "Neural", "Proneural"), 
                                          times = c(length(S1),length(S2), length(S3), length(S4)))))


names(dataset1) <- c('Drug sensitivity score', 'Subtype')


##### #####   Linsitinib_1510

dev.new()
ggboxplot(dataset1, x='Subtype', y='Drug sensitivity score', fill='Subtype', palette='jama',
          add='none'
)+
  geom_signif(comparisons = list(c("Mesenchymal", "Proneural"),c("Classical","Proneural"),c("Neural","Proneural")),
              y_position = c(180,190,200),
              tip_length = c(0),
              map_signif_level = F,
              test = wilcox.test
  )


##### #####  NVP-ADW742_1932

dev.new()
ggboxplot(dataset1, x='Subtype', y='Drug sensitivity score', fill='Subtype', palette='jama',
          add='none'
)+
  geom_signif(comparisons = list(c("Mesenchymal", "Proneural"),c("Classical","Proneural"),c("Neural","Proneural")),
              y_position = c(56,59,62),
              tip_length = c(0),
              map_signif_level = F,
              test = wilcox.test
  )


##### #####  GSK1904529A_1093

dev.new()
ggboxplot(dataset1, x='Subtype', y='Drug sensitivity score', fill='Subtype', palette='jama',
          add='none'
)+
  geom_signif(comparisons = list(c("Mesenchymal", "Proneural"),c("Classical","Proneural"),c("Neural","Proneural")),
              y_position = c(181,189,197),
              tip_length = c(0),
              map_signif_level = F,
              test = wilcox.test
  )


