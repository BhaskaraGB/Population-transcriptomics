library(ggplot2)
library(devtools)
library("BiocParallel")
library("vsn")
library(edgeR)
library(xlsx)

rm(list=ls())
register(MulticoreParam(workers=8))

#load("PH_Trans_01_13_20.RData")

info<-read.delim("/home/taslima/Tools/hmmer-3.2.1/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")
info<-info[,c(2,13,16)]
info<-info[!duplicated(info$locusName),]
dim(info)

####### Split 
### HAL_I vs. HAL_II
load("Split_H1_H2.RData")

sigGL_H_G_H1_H2<-GL_H_G_H1_H2 %>% filter( padj < (0.05/1))
sigGL_H_E_B_K<-GL_H_E_B_K %>%  filter( padj < (0.05/1))
sigGL_H_GE_B_K_H1_H2<-GL_H_GE_B_K_H1_H2 %>%  filter(padj <(0.05/1))

GGene<-merge(sigGL_H_G_H1_H2, info, by.x="row", by.y="locusName",all.y=F)
EGene<-merge(sigGL_H_E_B_K, info, by.x="row", by.y="locusName",all.y=F)
GEGene<-merge(sigGL_H_GE_B_K_H1_H2, info, by.x="row", by.y="locusName",all.y=F)

write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALI_vs_HALII.xlsx",x = GGene,sheetName = "GGene",append = F,row.names = F)
write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALI_vs_HALII.xlsx",x = EGene,sheetName = "EGene",append = T,row.names = F)
write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALI_vs_HALII.xlsx",x = GEGene,sheetName = "GEGene",append = T,row.names = F)

### HAL_I vs. FIL
rm(list=ls(pattern = "GL"))
rm(list=ls(pattern ="sig"))
load("Split_H1_F.RData")
ls()

sigGL_H1F_G_H1_F<-GL_H1F_G_H1_F %>% filter( padj < (0.05/1))
sigGL_H1F_E_B_K<-GL_H1F_E_B_K %>%  filter( padj < (0.05/1))
sigGL_H1F_GE_B_K_H1_F<-GL_H1F_GE_B_K_H1_F %>%  filter(padj <(0.05/1))

GGene<-merge(sigGL_H1F_G_H1_F, info, by.x="row", by.y="locusName",all.y=F)
EGene<-merge(sigGL_H1F_E_B_K, info, by.x="row", by.y="locusName",all.y=F)
GEGene<-merge(sigGL_H1F_GE_B_K_H1_F, info, by.x="row", by.y="locusName",all.y=F)

write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALI_vs_FIL.xlsx",x = GGene,sheetName = "GGene",append = F,row.names = F)
write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALI_vs_FIL.xlsx",x = EGene,sheetName = "EGene",append = T,row.names = F)
write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALI_vs_FIL.xlsx",x = GEGene,sheetName = "GEGene",append = T,row.names = F)

### HAL_II vs. FIL
rm(list=ls(pattern = "GL"))
rm(list=ls(pattern ="sig"))
load("Split_H2_F.RData")
ls()

sigGL_H2F_G_H2_F<-GL_H2F_G_H2_F %>% filter( padj < (0.05/1))
sigGL_H2F_E_B_K<-GL_H2F_E_B_K %>%  filter( padj < (0.05/1))
sigGL_H2F_GE_B_K_H2_F<-GL_H2F_GE_B_K_H2_F %>%  filter(padj <(0.05/1))

GGene<-merge(sigGL_H2F_G_H2_F, info, by.x="row", by.y="locusName",all.y=F)
EGene<-merge(sigGL_H2F_E_B_K, info, by.x="row", by.y="locusName",all.y=F)
GEGene<-merge(sigGL_H2F_GE_B_K_H2_F, info, by.x="row", by.y="locusName",all.y=F)

write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALII_vs_FIL.xlsx",x = GGene,sheetName = "GGene",append = F,row.names = F)
write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALII_vs_FIL.xlsx",x = EGene,sheetName = "EGene",append = T,row.names = F)
write.xlsx(file = "PH_Trans_3Grp_DEG_cont_HALII_vs_FIL.xlsx",x = GEGene,sheetName = "GEGene",append = T,row.names = F)


################ GO Enrichment
#info<-read.delim("/home/taslima/Tools/hmmer-3.2.1/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")

#Select only columns for geneID and GO
#info<-info[,c(2,10)]

#Get only genes that have GO annotation
#length(which(!is.na(info$GO)))

#info<-info[which(!is.na(info$GO)),]

#write.table(as.data.frame(info),file = "PH_HAL_v2_GO.tab",quote = F,row.names = F,col.names = F,sep = "\t")
#Install and load topGO
#source("https://bioconductor.org/biocLite.R")
#biocLite("topGO")
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")

library(topGO)
library(Rgraphviz)

#read customized annotation file
geneID2GO <- readMappings(file = "PH_HAL_v2_GO.tab")

## Test with random p-value
geneNames <- info$locusName

library(readxl)

## Read genelist from file
dat<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALI_vs_HALII.xlsx",sheetIndex = 1)

myInterestingGenes<-as.data.frame(dat[,1])
colnames(myInterestingGenes)<-"gene"
#myInterestingGenes <- targets$gene
geneList <- factor(as.integer(geneNames %in% myInterestingGenes$gene))
names(geneList) <- geneNames
str(geneList)

## BP
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=10)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
pvalFis <- score(resultFisher)
geneData(resultFisher)
allRes_BP <- GenTable(GOdata, classic = resultFisher, ranksOf = "classic",topNodes=as.vector(attr(resultFisher,"geneData"))[4],numChar=100)
allRes_BP$padj<-p.adjust(as.numeric(allRes_BP$classic),method = "fdr")
allRes_BP<-allRes_BP[which(allRes_BP$padj < as.vector(quantile(allRes_BP$padj,prob=0.1))),] # 1% FDR
allRes_BP$Attribute<-rep("BP",dim(allRes_BP)[1])

## MF

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=10)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
pvalFis <- score(resultFisher)
geneData(resultFisher)
allRes_MF <- GenTable(GOdata, classic = resultFisher, ranksOf = "classic",topNodes=as.vector(attr(resultFisher,"geneData"))[4],numChar=100)
allRes_MF$padj<-p.adjust(as.numeric(allRes_MF$classic),method = "fdr")
allRes_MF<-allRes_MF[which(allRes_MF$padj < as.vector(quantile(allRes_MF$padj,prob=0.1))),]
allRes_MF$Attribute<-rep("MF",dim(allRes_MF)[1])

##CC
GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=10)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes_CC <- GenTable(GOdata, classic = resultFisher, ranksOf = "classic",topNodes=as.vector(attr(resultFisher,"geneData"))[4],numChar=100)
allRes_CC$padj<-p.adjust(as.numeric(allRes_CC$classic),method = "fdr")
allRes_CC<-allRes_CC[which(allRes_MF$padj < as.vector(quantile(allRes_CC$padj,prob=0.1))),]
allRes_CC$Attribute<-rep("CC",dim(allRes_CC)[1])

write.xlsx(rbind(allRes_BP,allRes_CC,allRes_MF),file = "Test.xlsx",sheetName = "GO_Enrichment",row.names = F) ## U can also write in different sheet


