setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/PHtrans/")


library(ggplot2)
library(devtools)
library("BiocParallel")
library("vsn")
library(edgeR)
library(xlsx)
library(tidyverse)
source("Rscripts/RuntopGo.R")

rm(list=ls())

#Use multiple cores
register(MulticoreParam(workers=6))
# Set the memory for Java
options(java.parameters = "-Xmx6g" )

##### Read the annotation
info<-read.delim("DBs/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")
info<-info[,c(2,13,16)]
info<-info[!duplicated(info$locusName),]
dim(info)

########## Global GXE

load("Results/PH_Trans_Global.RData")

# Parse out by factor
sigGL_G <- GL_G %>% filter( padj < (0.05/1))
sigGL_E <- GL_E %>%  filter( padj < (0.05/1))
sigGL_GXE <- GL_GE %>%  filter(padj <(0.05/1))

# Four groups: Only G, Only E, G plus E (additive) and GxE (interactive)
onlyG<-setdiff(sigGL_G$row, unique(c(sigGL_E$row, sigGL_GXE$row)))
onlyE<-setdiff(sigGL_E$row,unique(c(sigGL_G$row, sigGL_GXE$row)))
GAndE<-setdiff( intersect(sigGL_E$row,sigGL_G$row), sigGL_GXE$row)
GXE<-sigGL_GXE$row

GGene<-merge(sigGL_G[which(sigGL_G$row %in% onlyG),], info, by.x="row", by.y="locusName",all.x=T, all.y=F)
EGene<-merge(sigGL_E[which(sigGL_E$row %in% onlyE),], info, by.x="row", by.y="locusName",all.x=T, all.y=F)
GAndEgene<-merge(sigGL_G[which(sigGL_G$row %in% GAndE),], info, by.x="row", by.y="locusName",all.x=T, all.y=F)
GXEGene<-merge(sigGL_GXE, info, by.x="row", by.y="locusName",all.y=F)

xlsx::write.xlsx(file = "Results/PH_Trans_4Grp_DEG_global.xlsx",x = GGene,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_4Grp_DEG_global.xlsx",x = EGene,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_4Grp_DEG_global.xlsx",x = GAndEgene,sheetName = "G+EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_4Grp_DEG_global.xlsx",x = GXEGene,sheetName = "GXEGene",append = T,row.names = F)

#### 3GRP

GeneG<-merge(sigGL_G, info, by.x="row", by.y="locusName",all.x=T, all.y=F)
GeneE<-merge(sigGL_E, info, by.x="row", by.y="locusName",all.x=T, all.y=F)
GeneGXE<-merge(sigGL_GXE, info, by.x="row", by.y="locusName",all.x=T, all.y=F)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_global.xlsx",x = GeneG,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_global.xlsx",x = GeneE,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_global.xlsx",x = GeneGXE,sheetName = "GXEGene",append = T,row.names = F)


########### Post hoc

####### H1vsH2
sigGL_H_G_H1_H2<- GL_G_H1_H2 %>% filter( padj < (0.05/1))
sigGL_H_E_B_K<-GL_H_E_B_K %>%  filter( padj < (0.05/1))
sigGL_H_GE_B_K_H1_H2<-GL_GE_B_K_H1_H2%>%  filter(padj <(0.05/1))

GGene<-merge(sigGL_H_G_H1_H2, info, by.x="row", by.y="locusName",all.y=F)
GEGene<-merge(sigGL_H_GE_B_K_H1_H2, info, by.x="row", by.y="locusName",all.y=F)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",x = GGene,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",x = GEGene,sheetName = "GEGene",append = T,row.names = F)


### HAL_I vs. FIL
sigGL_H1F_G_H1_F<-GL_G_H1_F %>% filter( padj < (0.05/1))
sigGL_H1F_GE_B_K_H1_F<-GL_GE_B_K_H1_F %>%  filter(padj <(0.05/1))

GGene<-merge(sigGL_H1F_G_H1_F, info, by.x="row", by.y="locusName",all.y=F)
GEGene<-merge(sigGL_H1F_GE_B_K_H1_F, info, by.x="row", by.y="locusName",all.y=F)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",x = GGene,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",x = GEGene,sheetName = "GEGene",append = T,row.names = F)

### HAL_II vs. FIL


sigGL_H2F_G_H2_F<-GL_G_H2_F %>% filter( padj < (0.05/1))
sigGL_H2F_GE_B_K_H2_F<-GL_GE_B_K_H2_F %>%  filter(padj <(0.05/1))

GGene<-merge(sigGL_H2F_G_H2_F, info, by.x="row", by.y="locusName",all.y=F)
GEGene<-merge(sigGL_H2F_GE_B_K_H2_F, info, by.x="row", by.y="locusName",all.y=F)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",x = GGene,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",x = GEGene,sheetName = "GEGene",append = T,row.names = F)


################ GO Enrichment
####### Read only the subset that is expressed in leaf for this experiment
geneID2GO <- "DBs/PH_HAL_v2_GO_subset.tab"
info<-read.delim("DBs/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")
head(info)
#Select only columns for geneID and GO
info<-info[,c(2,10,16)]
#Get only genes that have GO annotation
info<-info[which(!is.na(info$GO)),]

###### Requires "RunTopGO" function
GDEGH1F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",sheet = 1)
SigGOH1FG<-RunTopGO(geneID2GO,info,GDEGH1F$row)
EDEGH1F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",sheet = 2)
SigGOH1FE<-RunTopGO(geneID2GO,info,EDEGH1F$row)
GEDEGH1F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",sheet = 3)
SigGOH1FGE<-RunTopGO(geneID2GO,info,GEDEGH1F$row)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_FIL.xlsx",x = SigGOH1FG,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_FIL.xlsx",x = SigGOH1FE,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_FIL.xlsx",x = SigGOH1FGE,sheetName = "GEGene",append = T,row.names = F)

GDEGH2F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",sheet = 1)
SigGOH2FG<-RunTopGO(geneID2GO,info,GDEGH2F$row)
EDEGH2F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",sheet = 2)
SigGOH2FE<-RunTopGO(geneID2GO,info,EDEGH2F$row)
GEDEGH2F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",sheet = 3)
SigGOH2FGE<-RunTopGO(geneID2GO,info,GEDEGH2F$row)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALII_vs_FIL.xlsx",x = SigGOH2FG,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALII_vs_FIL.xlsx",x = SigGOH2FE,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALII_vs_FIL.xlsx",x = SigGOH2FGE,sheetName = "GEGene",append = T,row.names = F)

GDEGH1H2<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",sheet = 1)
SigGOH1H2G<-RunTopGO(geneID2GO,info,GDEGH1H2$row)
EDEGH1H2<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",sheet = 2)
SigGOH1H2E<-RunTopGO(geneID2GO,info,EDEGH1H2$row)
GEDEGH1H2<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",sheet = 3)
SigGOH1H2GE<-RunTopGO(geneID2GO,info,GEDEGH1H2$row)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_HALII.xlsx",x = SigGOH1H2G,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_HALII.xlsx",x = SigGOH1H2E,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_HALII.xlsx",x = SigGOH1H2GE,sheetName = "GEGene",append = T,row.names = F)



GDEGH1F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",sheet = 1)
SigGOH1FG<-RunTopGO(geneID2GO,info,GDEGH1F$row)
EDEGH1F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",sheet = 2)
SigGOH1FE<-RunTopGO(geneID2GO,info,EDEGH1F$row)
GEDEGH1F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_FIL.xlsx",sheet = 3)
SigGOH1FGE<-RunTopGO(geneID2GO,info,GEDEGH1F$row)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_FIL.xlsx",x = SigGOH1FG,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_FIL.xlsx",x = SigGOH1FE,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_FIL.xlsx",x = SigGOH1FGE,sheetName = "GEGene",append = T,row.names = F)

GDEGH2F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",sheet = 1)
SigGOH2FG<-RunTopGO(geneID2GO,info,GDEGH2F$row)
EDEGH2F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",sheet = 2)
SigGOH2FE<-RunTopGO(geneID2GO,info,EDEGH2F$row)
GEDEGH2F<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALII_vs_FIL.xlsx",sheet = 3)
SigGOH2FGE<-RunTopGO(geneID2GO,info,GEDEGH2F$row)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALII_vs_FIL.xlsx",x = SigGOH2FG,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALII_vs_FIL.xlsx",x = SigGOH2FE,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALII_vs_FIL.xlsx",x = SigGOH2FGE,sheetName = "GEGene",append = T,row.names = F)

GDEGH1H2<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",sheet = 1)
SigGOH1H2G<-RunTopGO(geneID2GO,info,GDEGH1H2$row)
EDEGH1H2<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",sheet = 2)
SigGOH1H2E<-RunTopGO(geneID2GO,info,EDEGH1H2$row)
GEDEGH1H2<-openxlsx::read.xlsx("Results/PH_Trans_3Grp_DEG_cont_Global_HALI_vs_HALII.xlsx",sheet = 3)
SigGOH1H2GE<-RunTopGO(geneID2GO,info,GEDEGH1H2$row)

xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_HALII.xlsx",x = SigGOH1H2G,sheetName = "GGene",append = F,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_HALII.xlsx",x = SigGOH1H2E,sheetName = "EGene",append = T,row.names = F)
xlsx::write.xlsx(file = "Results/PH_Trans_3Grp_GO_Global_HALI_vs_HALII.xlsx",x = SigGOH1H2GE,sheetName = "GEGene",append = T,row.names = F)


