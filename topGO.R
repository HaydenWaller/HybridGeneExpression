#topGO to identify enrichment

library(topGO)
library(Rgraphviz)
library(tidyverse)


# FAST vs SLOW DE

GOuniverse<-read_csv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/all.GOs.csv")
geneID2GO<-readMappings(file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/all.GOs.csv", sep = ",", IDsep = ";")
geneUniverse <-names(geneID2GO)
DE_genes<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/fast.slow.DE_GOs.csv", sep = ",", header = TRUE)
geneList<-factor(as.integer(geneUniverse %in% as.vector(DE_genes$Entry_Name)))
names(geneList)<-geneUniverse
myGOdata_BP<-new("topGOdata", description= "fastSlow", ontology = "BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata_MF<-new("topGOdata", description= "fastSlow", ontology = "MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
myGOdata_CC<-new("topGOdata", description= "fastSlow", ontology = "CC", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
resultFisher_BP <- runTest(myGOdata_BP, algorithm="weight01", statistic="fisher")
resultFisher_MF <- runTest(myGOdata_MF, algorithm="weight01", statistic="fisher")
resultFisher_CC <- runTest(myGOdata_CC, algorithm="weight01", statistic="fisher")

allRes_BP <- GenTable(myGOdata_BP, Fisher = resultFisher_BP, topNodes=150)
allRes_BP <- allRes_BP[which(allRes_BP$Significant > 0),]
allRes_BP$fdr_pvalue<-p.adjust(allRes_BP$Fisher,"fdr",nrow(allRes_BP))
allRes_MF <- GenTable(myGOdata_MF, Fisher = resultFisher_MF, topNodes=100)
allRes_MF <- allRes_MF[which(allRes_MF$Significant > 0),]
allRes_MF$fdr_pvalue<-p.adjust(allRes_MF$Fisher,"fdr",nrow(allRes_MF))
allRes_CC <- GenTable(myGOdata_CC, Fisher = resultFisher_CC, topNodes=100)
allRes_CC <- allRes_CC[which(allRes_CC$Significant > 0),]
allRes_CC$fdr_pvalue<-p.adjust(allRes_CC$Fisher,"fdr",nrow(allRes_CC))

write_tsv(allRes_MF, "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/fast.slow_MF.tsv")
write_tsv(allRes_BP, "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/fast.slow_BP.tsv")
write_tsv(allRes_CC, "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/fast.slow_CC.tsv")



# LK vs LP DE

GOuniverse<-read_csv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/all.GOs.csv")
geneID2GO<-readMappings(file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/all.GOs.csv", sep = ",", IDsep = ";")
geneUniverse <-names(geneID2GO)
DE_genes<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/lk.lp.DE_GOs.csv", sep = ",", header = TRUE)
geneList<-factor(as.integer(geneUniverse %in% as.vector(DE_genes$Entry_Name)))
names(geneList)<-geneUniverse
myGOdata_BP<-new("topGOdata", description= "lklp", ontology = "BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
myGOdata_MF<-new("topGOdata", description= "lklp", ontology = "MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
myGOdata_CC<-new("topGOdata", description= "lklp", ontology = "CC", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10) 
resultFisher_BP <- runTest(myGOdata_BP, algorithm="weight01", statistic="fisher")
resultFisher_MF <- runTest(myGOdata_MF, algorithm="weight01", statistic="fisher")
resultFisher_CC <- runTest(myGOdata_CC, algorithm="weight01", statistic="fisher")

allRes_BP <- GenTable(myGOdata_BP, Fisher = resultFisher_BP, topNodes=150)
allRes_BP <- allRes_BP[which(allRes_BP$Significant > 0),]
allRes_BP$fdr_pvalue<-p.adjust(allRes_BP$Fisher,"fdr",nrow(allRes_BP))
allRes_MF <- GenTable(myGOdata_MF, Fisher = resultFisher_MF, topNodes=100)
allRes_MF <- allRes_MF[which(allRes_MF$Significant > 0),]
allRes_MF$fdr_pvalue<-p.adjust(allRes_MF$Fisher,"fdr",nrow(allRes_MF))
allRes_CC <- GenTable(myGOdata_CC, Fisher = resultFisher_CC, topNodes=100)
allRes_CC <- allRes_CC[which(allRes_CC$Significant > 0),]
allRes_CC$fdr_pvalue<-p.adjust(allRes_CC$Fisher,"fdr",nrow(allRes_CC))

write_tsv(allRes_MF, "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/lk.lp_MF.tsv")
write_tsv(allRes_BP, "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/lk.lp_BP.tsv")
write_tsv(allRes_CC, "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/GO_Terms/lk.lp_CC.tsv")


