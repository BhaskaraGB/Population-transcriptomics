setwd("/home/taslima/data/JuengerLab/ALL_RNASeq/Bhashkar_RNASeq/V5")

rm(list=ls())

library(DESeq2)
library(ggplot2)
library(devtools)
library("BiocParallel")
library("vsn")
library(edgeR)
library(xlsx)

register(MulticoreParam(workers=24))

Design<-read.csv("PHTrans_Final_Meta.csv",row.names = 6)
dat<-read.table("PHTrans_Count_HALREF.tab",header = T,stringsAsFactors = F,check.names = F)

dat<-dat[,rownames(Design)]
dat<-dat[-which(as.vector(rowMeans(dat,na.rm = F)) <1),]  # Remove very low abundant gene from count matrix

dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = Design,
                              design = ~ 1  )

dds = estimateSizeFactors(dds )
sizeFactors(dds)

hist(as.vector(sizeFactors(dds)),breaks=50) # Check the distribution of library size

####### DEG Analysis

Design$SITE<-factor(Design$SITE,levels = c("KINGS","BFL"))
length(levels(Design$SITE))
Design$Population<-factor(Design$Population,levels = c("FIL","HAL_I","HAL_II"))
length(levels(Design$Population))

se<-SummarizedExperiment(assays = data.matrix(dat),
                         colData = DataFrame(Design))

flist = list(c("~ Population + SITE","~ SITE"),
             c("~ SITE + Population","~ Population"),
             c("~ Population + SITE + Population:SITE","~ Population + SITE"))

#geneIDs = rownames(dat)



##### G

x<-flist[[1]]
ddsG<- DESeqDataSet(se = se, design = as.formula(x[1]))
cat(paste(c("Start Running \"G Gene\" Analysis at",date(), "........\n")))
desG<-DESeq(ddsG,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desG)
GL_G_H1_H2<-results(desG,contrast = c("Population","HAL_I","HAL_II"),
        parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_G_H1_F<-results(desG,contrast = c("Population","HAL_I","FIL"),
                    parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)


GL_G_H2_F<-results(desG,contrast = c("Population","HAL_II","FIL"),
                   parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

#### E
x<-flist[[2]]
ddsE<- DESeqDataSet(se = se, design = as.formula(x[1]))
cat(paste(c("Start Running \"E Gene\" Analysis at",date(), "........\n")))
desE<-DESeq(ddsE,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desE)
GL_E_B_K<-results(desE,contrast = c("SITE","BFL","KINGS"),
                    parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)


#### GXE
x<-flist[[3]]
ddsGE<- DESeqDataSet(se = se, design = as.formula(x[1]))
cat(paste(c("Start Running \"GE Gene\" Analysis at",date(), "........\n")))
desGE<-DESeq(ddsGE,test="LRT", reduced=  as.formula(x[2]), parallel = T)

resultsNames(desGE)

GL_GE_B_K_H1_F<-results(desGE,name="PopulationHAL_I.SITEBFL",
                        parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_GE_B_K_H2_F<-results(desGE,name="PopulationHAL_II.SITEBFL",
                        parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_GE_B_K_H1_H2<-results(desGE,contrast = list(c("PopulationHAL_I.SITEBFL","PopulationHAL_I.SITEBFL")),
                        parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

save.image("PH_Trans_01_13_20.RData")

###### Split for contrasts

##  HAL_I vs HAL_II
DesignH<-Design[which(Design$Population %in% c("HAL_I","HAL_II")),]
dat_H<-dat[,rownames(DesignH)]

DesignH$Population<-factor(DesignH$Population,levels = c("HAL_II","HAL_I"))
length(levels(DesignH$Population))
DesignH$SITE<-factor(DesignH$SITE,levels = c("KINGS","BFL"))
length(levels(DesignH$SITE))

seH<-SummarizedExperiment(assays = data.matrix(dat_H),
                         colData = DataFrame(DesignH))

flist = list(c("~ Population + SITE","~ SITE"),
             c("~ SITE + Population","~ Population"),
             c("~ Population + SITE + Population:SITE","~ Population + SITE"))

##### G

x<-flist[[1]]
ddsH_G<- DESeqDataSet(se = seH, design = as.formula(x[1]))
cat(paste(c("Start Running split for HAL_I vs. HAL_II \"G Gene\" Analysis at",date(), "........\n")))
desH_G<-DESeq(ddsH_G,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desH_G)
GL_H_G_H1_H2<-results(desH_G,contrast = c("Population","HAL_I","HAL_II"),
                    parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

##### E
x<-flist[[2]]
ddsH_E<- DESeqDataSet(se = seH, design = as.formula(x[1]))
cat(paste(c("Start Running split for HAL_I vs. HAL_II \"E Gene\" Analysis at",date(), "........\n")))
desH_E<-DESeq(ddsH_E,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desH_E)
GL_H_E_B_K<-results(desH_E,contrast = c("SITE","BFL","KINGS"),
                  parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

#### GXE
x<-flist[[3]]
ddsH_GE<- DESeqDataSet(se = seH, design = as.formula(x[1]))
cat(paste(c("Start Running split HAL_I vs. HAL_II \"GE Gene\" Analysis at",date(), "........\n")))
desH_GE<-DESeq(ddsH_GE,test="LRT", reduced=  as.formula(x[2]), parallel = T)

resultsNames(desH_GE)

GL_H_GE_B_K_H1_H2<-results(desH_GE,name = "PopulationHAL_I.SITEBFL",
                         parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)


save(GL_H_G_H1_H2,GL_H_E_B_K,GL_H_GE_B_K_H1_H2,file="Split_H1_H2.RData")

####################################
##  HAL_I vs FIL
DesignH1F<-Design[which(Design$Population %in% c("HAL_I","FIL")),]
dat_H1F<-dat[,rownames(DesignH1F)]

DesignH1F$Population<-factor(DesignH1F$Population,levels = c("FIL","HAL_I"))
length(levels(DesignH1F$Population))
DesignH1F$SITE<-factor(DesignH1F$SITE,levels = c("KINGS","BFL"))
length(levels(DesignH1F$SITE))

seH1F<-SummarizedExperiment(assays = data.matrix(dat_H1F),
                          colData = DataFrame(DesignH1F))

flist = list(c("~ Population + SITE","~ SITE"),
             c("~ SITE + Population","~ Population"),
             c("~ Population + SITE + Population:SITE","~ Population + SITE"))

##### G

x<-flist[[1]]
ddsH1F_G<- DESeqDataSet(se = seH1F, design = as.formula(x[1]))
cat(paste(c("Start Running split for HAL_I vs. FIL \"G Gene\" Analysis at",date(), "........\n")))
desH1F_G<-DESeq(ddsH1F_G,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desH1F_G)
GL_H1F_G_H1_F<-results(desH1F_G,contrast = c("Population","HAL_I","FIL"),
                      parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

##### E
x<-flist[[2]]
ddsH1F_E<- DESeqDataSet(se = seH1F, design = as.formula(x[1]))
cat(paste(c("Start Running split for HAL_I vs. FIL \"E Gene\" Analysis at",date(), "........\n")))
desH1F_E<-DESeq(ddsH1F_E,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desH1F_E)
GL_H1F_E_B_K<-results(desH1F_E,contrast = c("SITE","BFL","KINGS"),
                    parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

#### GXE
x<-flist[[3]]
ddsH1F_GE<- DESeqDataSet(se = seH1F, design = as.formula(x[1]))
cat(paste(c("Start Running split HAL_I vs. FIL \"GE Gene\" Analysis at",date(), "........\n")))
desH1F_GE<-DESeq(ddsH1F_GE,test="LRT", reduced=  as.formula(x[2]), parallel = T)

resultsNames(desH1F_GE)

GL_H1F_GE_B_K_H1_F<-results(desH1F_GE,name = "PopulationHAL_I.SITEBFL",
                           parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

save(GL_H1F_G_H1_F,GL_H1F_E_B_K,GL_H1F_GE_B_K_H1_F,file="Split_H1_F.RData")

####################################
##  HAL_II vs FIL
DesignH2F<-Design[which(Design$Population %in% c("HAL_II","FIL")),]
dat_H2F<-dat[,rownames(DesignH2F)]

DesignH2F$Population<-factor(DesignH2F$Population,levels = c("FIL","HAL_II"))
length(levels(DesignH2F$Population))
DesignH2F$SITE<-factor(DesignH2F$SITE,levels = c("KINGS","BFL"))
length(levels(DesignH2F$SITE))

seH2F<-SummarizedExperiment(assays = data.matrix(dat_H2F),
                            colData = DataFrame(DesignH2F))

flist = list(c("~ Population + SITE","~ SITE"),
             c("~ SITE + Population","~ Population"),
             c("~ Population + SITE + Population:SITE","~ Population + SITE"))

##### G

x<-flist[[1]]
ddsH2F_G<- DESeqDataSet(se = seH2F, design = as.formula(x[1]))
cat(paste(c("Start Running split for HAL_II vs. FIL \"G Gene\" Analysis at",date(), "........\n")))
desH2F_G<-DESeq(ddsH2F_G,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desH2F_G)
GL_H2F_G_H2_F<-results(desH2F_G,contrast = c("Population","HAL_II","FIL"),
                       parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

##### E
x<-flist[[2]]
ddsH2F_E<- DESeqDataSet(se = seH2F, design = as.formula(x[1]))
cat(paste(c("Start Running split for HAL_II vs. FIL \"E Gene\" Analysis at",date(), "........\n")))
desH2F_E<-DESeq(ddsH2F_E,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desH2F_E)
GL_H2F_E_B_K<-results(desH2F_E,contrast = c("SITE","BFL","KINGS"),
                      parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

#### GXE
x<-flist[[3]]
ddsH2F_GE<- DESeqDataSet(se = seH2F, design = as.formula(x[1]))
cat(paste(c("Start Running split HAL_II vs. FIL \"GE Gene\" Analysis at",date(), "........\n")))
desH2F_GE<-DESeq(ddsH2F_GE,test="LRT", reduced=  as.formula(x[2]), parallel = T)

resultsNames(desH2F_GE)

GL_H2F_GE_B_K_H2_F<-results(desH2F_GE,name = "PopulationHAL_II.SITEBFL",
                            parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

save(GL_H2F_G_H2_F,GL_H2F_E_B_K,GL_H2F_GE_B_K_H2_F,file="Split_H2_F.RData")
