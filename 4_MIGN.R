

install.packages("survival")
library(survival)

install.packages("CPE")
library(CPE)

library(ggpubr)
library(dplyr)

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


TAM_brown=read.csv("TAMgeneSet_2.csv")
TAM_brown=as.character(TAM_brown[,2])

TAM_turquoise=read.csv("TAMgeneSet_4.csv")
TAM_turquoise=as.character(TAM_turquoise[,2])

TC_blue=read.csv("TCgeneSet_1.csv")
TC_blue=as.character(TC_blue[,2])

TC_brown=read.csv("TCgeneSet_2.csv")
TC_brown=as.character(TC_brown[,2])

TC_green=read.csv("TCgeneSet_3.csv")
TC_green=as.character(TC_green[,2])

TC_turquoise=read.csv("TCgeneSet_5.csv")
TC_turquoise=as.character(TC_turquoise[,2])

TC_yellow=read.csv("TCgeneSet_6.csv")
TC_yellow=as.character(TC_yellow[,2])


Marker_genes=unique(c(TAM_brown, TAM_turquoise, TC_blue, TC_green, TC_turquoise))  
Marker_genes = Reduce(intersect, list(Marker_genes, rownames(Gene_CGGALGG), rownames(Gene_LGG)))  # 995


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


Gene_LGG=Gene_LGG[Marker_genes,]  # 881  525

datExpr0=matrix(as.numeric(Gene_LGG),dim(Gene_LGG))
datExpr0=t(datExpr0)
rownames(datExpr0)=colnames(Gene_LGG)
colnames(datExpr0)=Marker_genes

datTraits=Survival_LGG


###
Gene_CGGALGG=Gene_CGGALGG[Marker_genes,]  # 259  172

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
                                                    minClusterSize = 40))

# Now we plot the results
dev.new()
sizeGrWindow(8,6)
plotDendroAndColors(hierTOM, colorDynamicHybridTOM, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

table(colorDynamicHybridTOM)

#colorDynamicHybridTOM (blue, turquoise)
#blue     brown      grey turquoise 
#219        85         6       240 

#colorDynamicHybridTOM (blue, turquoise, brown)
#blue     brown      grey turquoise    yellow 
#206       191         8       237        78 

#colorDynamicHybridTOM (blue, turquoise, green)
#blue     brown      grey turquoise 
#220       149        14       226 

#colorDynamicHybridTOM (blue, turquoise, yellow)
#blue     brown      grey turquoise    yellow 
#218       109         6       232        79


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


colorh1= colorDynamicHybridTOM
datME=moduleEigengenes(datExpr,colorh1)$eigengenes


table(colorh1)
#blue     brown      grey turquoise 
#220       149        14       226 

geneNames=colnames(datExpr)

geneSet_1=geneNames[colorh1=="blue"]
geneSet_2=geneNames[colorh1=="brown"]
geneSet_3=geneNames[colorh1=="grey"]
geneSet_4=geneNames[colorh1=="turquoise"]

write.csv(geneSet_1, file = "geneSet_1.csv")
write.csv(geneSet_2, file = "geneSet_2.csv")
write.csv(geneSet_3, file = "geneSet_3.csv")
write.csv(geneSet_4, file = "geneSet_4.csv")


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

#blue     brown      grey  turquoise 
#220       149        14       226 
#0.673    0.623      0.591 0.656


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


dataset1 <- data.frame(value = c(GS1,GS2,GS3,GS4), 
                       group = factor(rep(c("Blue", "Brown","Grey","Turquoise"), 
                                          times = c(length(GS1),length(GS2),length(GS3),length(GS4)))))


names(dataset1) <- c('K-index for overall survival', 'Gene Module')

dev.new()
ggboxplot(dataset1, x='Gene Module', y='K-index for overall survival', fill='Gene Module', 
          add='point')+
  scale_fill_manual(values = c("blue","brown","grey","turquoise"))+
  stat_compare_means(method = "anova", label.y = 0.74)+
  stat_compare_means(label = "p.signif", method = "t.test", ref.group=c("Blue"), label.y = 0.73)+
  stat_compare_means(label = "p.signif", method = "t.test", ref.group=c("Brown"), label.y = 0.72)+
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
colorlevels=c('blue','brown','grey','turquoise')
sizeGrWindow(12,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]]; 
  restrict1 = (colorh1==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1], 
                     GeneSignificance[restrict1], col=colorh1[restrict1],
                     main=whichmodule,  ylim = c(0.5,0.7),
                     xlab = "Connectivity", ylab = "K-index", abline = TRUE,
                     cex=2.2, pch=20)
}


#=====================================================================================
#
#  Code chunk 18
#
#=====================================================================================


### ### blue
Alldegrees_blue=Alldegrees1[colorh1=="blue",]
kTotoal_blue=Alldegrees_blue$kWithin
names(kTotoal_blue)=rownames(Alldegrees_blue)

FilterGenes= GS1>0.6 & kTotoal_blue > quantile(kTotoal_blue,0.9)
table(FilterGenes)
selected_blue=names(GS1)[FilterGenes]  # 21


### ### brown
Alldegrees_brown=Alldegrees1[colorh1=="brown",]
kTotoal_brown=Alldegrees_brown$kWithin
names(kTotoal_brown)=rownames(Alldegrees_brown)

FilterGenes= GS2>0.6 & kTotoal_brown > quantile(kTotoal_brown,0.9)
table(FilterGenes)
selected_brown=names(GS2)[FilterGenes]  # 6


### ### turquoise
Alldegrees_tur=Alldegrees1[colorh1=="turquoise",]
kTotoal_tur=Alldegrees_tur$kWithin
names(kTotoal_tur)=rownames(Alldegrees_tur)

FilterGenes= GS4>0.6 & kTotoal_tur > quantile(kTotoal_tur,0.9)
table(FilterGenes)
selected_tur=names(GS4)[FilterGenes]  # 23


hub_gene=c(selected_blue, selected_brown, selected_tur)
saveRDS(hub_gene, "hub_gene.rds")


#=====================================================================================
#
#  Code chunk 19
#
#=====================================================================================


LigRec = read.table("database/LigRec.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)


Marker_TAM=read.csv("humanTAM.csv")
Marker_TAM=as.character(Marker_TAM[,2])

Marker_TC=read.csv("humanTC.csv")
Marker_TC=as.character(Marker_TC[,2])

Ligand_blue=intersect(unique(LigRec$Ligand),geneSet_1)       #9
Receptor_blue=intersect(unique(LigRec$Receptor),geneSet_1)   #23

Ligand_brown=intersect(unique(LigRec$Ligand),geneSet_2)      #19
Receptor_brown=intersect(unique(LigRec$Receptor),geneSet_2)  #23


Ligand = unique(c(intersect(Ligand_blue,Marker_TAM), 
                        intersect(Marker_TAM,Ligand_brown)))  #13

Receptor = unique(c(intersect(Receptor_blue,Marker_TC),
                        intersect(Marker_TC,Receptor_brown))) #27


saveRDS(Ligand, "Ligand.rds")

saveRDS(Receptor, "Receptor.rds")





