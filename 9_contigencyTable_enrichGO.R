

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


TAM_brown=read.csv("TAMgeneSet_2.csv")
TAM_brown=as.character(TAM_brown[,2])

TAM_turquoise=read.csv("TAMgeneSet_4.csv")
TAM_turquoise=as.character(TAM_turquoise[,2])

TC_blue=read.csv("TCgeneSet_1.csv")
TC_blue=as.character(TC_blue[,2])

TC_green=read.csv("TCgeneSet_3.csv")
TC_green=as.character(TC_green[,2])

TC_turquoise=read.csv("TCgeneSet_5.csv")
TC_turquoise=as.character(TC_turquoise[,2])

MIGN_blue = read.csv("geneSet_1.csv")
MIGN_blue = as.character(MIGN_blue[,2])

MIGN_brown = read.csv("geneSet_2.csv")
MIGN_brown = as.character(MIGN_brown[,2])

MIGN_turquoise = read.csv("geneSet_4.csv")
MIGN_turquoise = as.character(MIGN_turquoise[,2])


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


TAM_gene = union(TAM_brown, TAM_turquoise)

#MIGN = list(MIGN_blue, MIGN_brown=MIGN_brown, MIGN_turquoise=MIGN_turquoise)
#TAM = list(TAM_brown=TAM_brown, TAM_turquoise=TAM_turquoise)

MIGN = list(MIGN_blue, MIGN_brown, MIGN_turquoise)
TAM = list(TAM_brown, TAM_turquoise)

MIGN_color = c("blue", "brown", "turquoise")
TAM_color = c("brown", "turquoise")



# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = 3, ncol = 2);
# Execute all pairwaise comparisons
CountTbl = matrix(0, nrow = 3, ncol = 2);


for (lmod in 1:2)
  for (smod in 1:3)
  {
    TAMMmember = (TAM_gene %in% MIGN[[smod]]);
    TAMHmember = (TAM_gene %in% TAM[[lmod]]);
    pTable[smod,lmod] = -log10(fisher.test(TAMMmember, TAMHmember, alternative = "greater")$p.value);
    CountTbl[smod,lmod] = length(intersect(MIGN[[smod]],TAM[[lmod]]))
  }

# Create text matrix
textMatrix = paste(CountTbl, "\n(",
                   signif(pTable, 3), ")", sep = "");
# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
TAMMModTotals = c(220,149,226)
TAMHModTotals = c(38,117)
# Actual plotting
dev.new()
sizeGrWindow(10,7 );
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
# jpeg(file = "Multicellular network/Network construction Proneural/Figure/Overlap_TAM.jpeg", wi = 800, he = 600);
#cairo_ps("Multicellular network/Network construction Proneural/Figure/For thesis/Overlap_TAM.eps",width = 12,height=9)
par(mfrow=c(1,1));
par(cex = 1.1);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", TAM_color),
               yLabels = paste(" ", MIGN_color),
               colorLabels = TRUE,
               xSymbols = paste("TAM ", TAM_color, ": ", TAMHModTotals, sep=""),
               ySymbols = paste("MIGN ", MIGN_color, ": ", TAMMModTotals, sep=""),
               textMatrix = textMatrix, # CountTbl
               colors = blueWhiteRed(100)[50:100],
               cex.text = 1.1, cex.lab = 1.1, setStdMargins = FALSE)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


TC_gene = c(TC_blue, TC_green, TC_turquoise)

#MIGN = list(MIGN_blue, MIGN_brown=MIGN_brown, MIGN_turquoise=MIGN_turquoise)
#TAM = list(TAM_brown=TAM_brown, TAM_turquoise=TAM_turquoise)

MIGN = list(MIGN_blue, MIGN_brown, MIGN_turquoise)
TC = list(TC_blue, TC_green, TC_turquoise)

MIGN_color = c("blue", "brown", "turquoise")
TC_color = c("blue", "green", "turquoise")



# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = 3, ncol = 3);
# Execute all pairwaise comparisons
CountTbl = matrix(0, nrow = 3, ncol = 3);


for (lmod in 1:3)
  for (smod in 1:3)
  {
    TAMMmember = (TC_gene %in% MIGN[[smod]]);
    TAMHmember = (TC_gene %in% TC[[lmod]]);
    pTable[smod,lmod] = -log10(fisher.test(TAMMmember, TAMHmember, alternative = "greater")$p.value);
    CountTbl[smod,lmod] = length(intersect(MIGN[[smod]],TC[[lmod]]))
  }

# Create text matrix
textMatrix = paste(CountTbl, "\n(",
                   signif(pTable, 3), ")", sep = "");
# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
TAMMModTotals = c(220,149,226)
TAMHModTotals = c(182,65,223)
# Actual plotting
dev.new()
sizeGrWindow(10,7 );
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
# jpeg(file = "Multicellular network/Network construction Proneural/Figure/Overlap_TAM.jpeg", wi = 800, he = 600);
#cairo_ps("Multicellular network/Network construction Proneural/Figure/For thesis/Overlap_TAM.eps",width = 12,height=9)
par(mfrow=c(1,1));
par(cex = 1.1);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", TC_color),
               yLabels = paste(" ", MIGN_color),
               colorLabels = TRUE,
               xSymbols = paste("TC ", TC_color, ": ", TAMHModTotals, sep=""),
               ySymbols = paste("MIGN ", MIGN_color, ": ", TAMMModTotals, sep=""),
               textMatrix = textMatrix, # CountTbl
               colors = blueWhiteRed(100)[50:100],
               cex.text = 1.1, cex.lab = 1.1, setStdMargins = FALSE)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


source("https://bioconductor.org/biocLite.R")

BiocManager::install("clusterProfiler")  #用来做富集分析
BiocManager::install("topGO")  #画GO图用的
BiocManager::install("Rgraphviz")
BiocManager::install("pathview") #看KEGG pathway的
BiocManager::install("org.Hs.eg.db") #这个包里存有人的注释文件
BiocManager::install("DOSE") 

# 载入包dian
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(DOSE)


##### #####


Marker_genes=MIGN_blue

Marker_genes=MIGN_brown


DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = Marker_genes,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

DEG.entrez_id = na.omit(DEG.entrez_id)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


##### GO enrichment (biological process) ##### (1)
erich.go.BP = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dev.new()
dotplot(erich.go.BP)

#pdf(file="./enrich.go.bp.tree.pdf",width = 10,height = 15)
#plotGOgraph(erich.go.BP)
#dev.off()



##### GO enrichment (Molecular function) ##### (2)***
erich.go.MF = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dev.new()
dotplot(erich.go.MF)



##### GO enrichment (Cellular component) ##### (3)
erich.go.CC = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

#barplot(erich.go.CC)

dev.new()
dotplot(erich.go.CC)



##### KEGG enrichment (biological process) ##### (4)
enrich.KEGG <- enrichKEGG(gene = DEG.entrez_id, organism = "hsa", 
                          keyType = "kegg", pvalueCutoff = 0.5, 
                          qvalueCutoff = 0.5)

dev.new()
dotplot(enrich.KEGG)

kegg=data.frame(enrich.KEGG)



##### DO enrichment ##### (5)
erich.do = enrichDO(gene = DEG.entrez_id,ont = "DO",pvalueCutoff = 0.5,qvalueCutoff = 0.5)

dev.new()
dotplot(erich.do)

signature_do=data.frame(erich.do)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


library(tidyverse)
library(org.Hs.eg.db)


GO_bp=data.frame(erich.go.BP)


### gliogenesis
DNA_geneID <- get('GO:0042063', org.Hs.egGO2ALLEGS) 
head(DNA_geneID)
length(DNA_geneID)   # 325

DNA_geneSYMBOL <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 

intersect(DNA_geneSYMBOL, MIGN_blue)
length(intersect(DNA_geneSYMBOL, MIGN_blue))  



### glial cell differentiation
DNA_geneID <- get('GO:0010001', org.Hs.egGO2ALLEGS) 
head(DNA_geneID)
length(DNA_geneID)   # 231

DNA_geneSYMBOL <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 



### regulation of nervous system development
DNA_geneID <- get('GO:0051960', org.Hs.egGO2ALLEGS) 
head(DNA_geneID)
length(DNA_geneID)   # 483

DNA_geneSYMBOL <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 



### positive regulation of tumor necrosis factor production
DNA_geneID <- get('GO:0032760', org.Hs.egGO2ALLEGS) 
head(DNA_geneID)
length(DNA_geneID)   # 96

DNA_geneSYMBOL <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 



### positive regulation of tumor necrosis factor superfamily cytokine production
DNA_geneID <- get('GO:1903557', org.Hs.egGO2ALLEGS) 
head(DNA_geneID)
length(DNA_geneID)   # 100

DNA_geneSYMBOL <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 



### glial cell development
DNA_geneID <- get('GO:0021782', org.Hs.egGO2ALLEGS) 
head(DNA_geneID)
length(DNA_geneID)   # 120

DNA_geneSYMBOL <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 



### glioblastoma multiforme
DNA_geneID <- get('DOID:3068', org.Hs.egDO2ALLEGS) 
head(DNA_geneID)
length(DNA_geneID)   # 120

DNA_geneSYMBOL <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 




