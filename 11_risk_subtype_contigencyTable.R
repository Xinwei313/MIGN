

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


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


groups = readRDS("groupsCGGA.rds")


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


lowRisk = sample[groups==1]
highRisk = sample[groups==2]


classical = sample[Survival_CGGA$Subtype=="Classical"]
mesenchymal = sample[Survival_CGGA$Subtype=="Mesenchymal"]
neural = sample[Survival_CGGA$Subtype=="Neural"]
proneural = sample[Survival_CGGA$Subtype=="Proneural"]


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sampleCGGA = sample

riskGroup = list(lowRisk, highRisk)
subtype = list(classical, mesenchymal, neural, proneural)

risk_color = c("blue", "brown")
subtype_color = c("brown", "turquoise","yellow","green")

subtype_group = c('Classical', 'Mesenchymal', 'Neural', 'Proneural')
risk_group = c('Low-risk group', 'High-risk group')

# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = 2, ncol = 4);
# Execute all pairwaise comparisons
CountTbl = matrix(0, nrow = 2, ncol = 4);


for (lmod in 1:4)
  for (smod in 1:2)
  {
    TAMMmember = (sampleCGGA %in% riskGroup[[smod]]);
    TAMHmember = (sampleCGGA %in% subtype[[lmod]]);
    pTable[smod,lmod] = -log10(fisher.test(TAMMmember, TAMHmember, alternative = "greater")$p.value);
    CountTbl[smod,lmod] = length(intersect(riskGroup[[smod]],subtype[[lmod]]))
  }

# Create text matrix
textMatrix = paste(CountTbl, "\n(",
                   signif(pTable, 3), ")", sep = "");
# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
TAMMModTotals = c(86,86)
TAMHModTotals = c(23,15,65,69)
# Actual plotting
dev.new()
sizeGrWindow(10,7);
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
# jpeg(file = "Multicellular network/Network construction Proneural/Figure/Overlap_TAM.jpeg", wi = 800, he = 600);
#cairo_ps("Multicellular network/Network construction Proneural/Figure/For thesis/Overlap_TAM.eps",width = 12,height=9)
par(mfrow=c(1,1));
par(cex = 1.1);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", subtype_color),
               yLabels = paste(" ", risk_color),
               colorLabels = TRUE,
               xSymbols = paste(subtype_group, ": ", TAMHModTotals, sep=""),
               ySymbols = paste(risk_group, ": ", TAMMModTotals, sep=""),
               textMatrix = textMatrix, # CountTbl
               colors = blueWhiteRed(100)[50:100],
               cex.text = 1.1, cex.lab = 1.1, setStdMargins = FALSE)

