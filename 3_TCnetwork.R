

install.packages("survival")
library(survival)

install.packages("CPE")
library(CPE)

library(ggpubr)


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


Marker_TC=read.csv("humanTC.csv")
Marker_TC=as.character(Marker_TC[,2])

Marker_genes = Reduce(intersect, list(Marker_TC,rownames(Gene_CGGALGG),rownames(Gene_LGG)))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


Gene_LGG=Gene_LGG[Marker_genes,]  # 784  525

datExpr0=matrix(as.numeric(Gene_LGG),dim(Gene_LGG))
datExpr0=t(datExpr0)
rownames(datExpr0)=colnames(Gene_LGG)
colnames(datExpr0)=Marker_genes

datTraits=Survival_LGG


###
Gene_CGGALGG=Gene_CGGALGG[Marker_genes,]  # 784  172

datExpr0=matrix(as.numeric(Gene_CGGALGG),dim(Gene_CGGALGG))
datExpr0=t(datExpr0)
rownames(datExpr0)=colnames(Gene_CGGALGG)
colnames(datExpr0)=Marker_genes

datTraits=Survival_CGGALGG


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GO.db", "preprocessCore", "impute")) 


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory.  On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load WGCNA package
library(WGCNA)
# Load additional necessary packages
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{ # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
dev.new()
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 60, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

datTraits=datTraits[rownames(datExpr),]

rownames(datExpr0)[clust==0]
#"TCGA-P5-A5F6-01"


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


dev.new()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^5
# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
#k=softConnectivity(datE=datExpr,power=5) 
# Plot a histogram of k and a scale free topology plot
dev.new()
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


dissTOM=TOMdist(ADJ1)

hierTOM = hclust(as.dist(dissTOM),method="average");

colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM = dissTOM, 
                                                    deepSplit = 3, pamRespectsDendro = FALSE,
                                                    minClusterSize = 30))

# Now we plot the results
dev.new()
sizeGrWindow(8,6)
plotDendroAndColors(hierTOM, colorDynamicHybridTOM, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

table(colorDynamicHybridTOM)

#colorDynamicHybridTOM (beta=5, deepSplit = 3, minClusterSize = 30)
#blue     brown     green      grey turquoise    yellow 
#182       172        65        48       223        94 


#=====================================================================================
#
#  Code chunk 12 (heatmap)
#
#=====================================================================================


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 5);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^5;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
dev.new()
sizeGrWindow(9,9)
#TOMplot(plotTOM, hierTOM, colorDynamicHybridTOM, main = "Network heatmap plot")
TOMplot(plotTOM, hierTOM, colorDynamicHybridTOM, main = "Network heatmap plot", col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))


#=====================================================================================
#
#  Code chunk 13
#
#=====================================================================================


colorh1 = colorDynamicHybridTOM
datME=moduleEigengenes(datExpr,colorh1)$eigengenes


table(colorh1)
#blue     brown     green      grey turquoise    yellow 
#182       172        65        48       223        94 

geneNames=colnames(datExpr)

geneSet_1=geneNames[colorh1=="blue"]
geneSet_2=geneNames[colorh1=="brown"]
geneSet_3=geneNames[colorh1=="green"]
geneSet_4=geneNames[colorh1=="grey"]
geneSet_5=geneNames[colorh1=="turquoise"]
geneSet_6=geneNames[colorh1=="yellow"]

write.csv(geneSet_1, file = "TCgeneSet_1.csv")
write.csv(geneSet_2, file = "TCgeneSet_2.csv")
write.csv(geneSet_3, file = "TCgeneSet_3.csv")
write.csv(geneSet_4, file = "TCgeneSet_4.csv")
write.csv(geneSet_5, file = "TCgeneSet_5.csv")
write.csv(geneSet_6, file = "TCgeneSet_6.csv")


#=====================================================================================
#
#  Code chunk 14
#
#=====================================================================================


MS=matrix(0,1,dim(datME)[2])

for (i in 1:dim(datME)[2])
{
  Gene_marker=as.matrix(datME)[,i]
  
  Survival=as.matrix(datTraits)
  
  status=as.numeric(Survival[,2])
  time=as.numeric(Survival[,3])
  
  fit <- coxph(Surv(time, status) ~ Gene_marker, method="breslow")
  
  #C-index
  #p_7[i]=summary(fit)$concordance[1]
  
  #K-index
  kindex=phcpe(coxfit=fit, CPE.SE = F)
  MS[i]=kindex$CPE
}
round(MS,3)


#blue     brown     green      grey turquoise    yellow 
#182       172        65        48       223        94 
#0.669    0.55      0.599     0.601  0.656       0.552


#=====================================================================================
#
#  Code chunk 15
#
#=====================================================================================


p=matrix(0,1,length(geneNames))

for (i in 1:length(geneNames))
{
  Marker_genes=c(geneNames[i])  
  
  Gene_marker=datExpr[,t(Marker_genes)]
  
  Survival=as.matrix(datTraits)
  
  status=as.numeric(Survival[,2])
  time=as.numeric(Survival[,3])
  
  fit <- coxph(Surv(time, status) ~ Gene_marker, method="breslow")
  
  #C-index
  #p[i]=summary(fit)$concordance[1]
  
  #K-index
  kindex=phcpe(coxfit=fit, CPE.SE = F)
  p[i]=kindex$CPE
}

GeneSignificance=as.double(p)
names(GeneSignificance)=colnames(datExpr)
GeneSignificance[1:5]


#=====================================================================================
#
#  Code chunk 16 (boxplot)
#
#=====================================================================================


GS1=GeneSignificance[geneSet_1]
GS2=GeneSignificance[geneSet_2]
GS3=GeneSignificance[geneSet_3]
GS4=GeneSignificance[geneSet_4]
GS5=GeneSignificance[geneSet_5]
GS6=GeneSignificance[geneSet_6]


dataset1 <- data.frame(value = c(GS1,GS2,GS3,GS4,GS5,GS6), 
                       group = factor(rep(c("Blue", "Brown","Green","Grey","Turquoise","Yellow"), 
                                          times = c(length(GS1),length(GS2),length(GS3),length(GS4),length(GS5),length(GS6)))))


names(dataset1) <- c('K-index for overall survival', 'Gene Module')

dev.new()
ggboxplot(dataset1, x='Gene Module', y='K-index for overall survival', fill='Gene Module', 
          add='point')+
  scale_fill_manual(values = c("blue","brown","green","grey","turquoise","yellow"))+
  stat_compare_means(method = "anova", label.y = 0.74)+
  stat_compare_means(label = "p.signif", method = "t.test", ref.group=c("Blue"), label.y = 0.73)+
  stat_compare_means(label = "p.signif", method = "t.test", ref.group=c("Green"), label.y = 0.72)+
  stat_compare_means(label = "p.signif", method = "t.test", ref.group=c("Turquoise"), label.y = 0.71)


#=====================================================================================
#
#  Code chunk 17 (scatter plot)
#
#=====================================================================================


ADJ1=abs(cor(datExpr,use="p"))^5
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)


dev.new()
#colorlevels=c("blue","turquoise")
colorlevels=c('blue','brown','green','grey','turquoise','yellow')
sizeGrWindow(12,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]]; 
  restrict1 = (colorh1==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1], 
                     GeneSignificance[restrict1], col=colorh1[restrict1],
                     ylim=c(0.5,0.7),
                     xlab = "Connectivity", ylab = "K-index", abline = TRUE,
                     cex=2.2, pch=20)
}

