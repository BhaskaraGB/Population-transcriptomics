setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/PHtrans/")

rm(list=ls())

library(DESeq2)
library(ggplot2)
library(devtools)
library("BiocParallel")
library("vsn")
library(edgeR)
library(xlsx)

register(MulticoreParam(workers=6))


# Read Design file
Design<-read.csv("Data/PHTrans_Final_Meta_548lib_FullStat.csv",row.names = 1)
#REad count matrix
dat<-read.table("Data/PHTrans_Count_HALREF_ALL.tab",header = T,stringsAsFactors = F,check.names = F)

# Filter low counts
dat<-dat[,rownames(Design)]
dat<-dat[-which(as.vector(rowMeans(dat,na.rm = F)) <1),]  # Remove very low abundant gene from count matrix

########### Global Analysis
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

GL_G<-results(desG,test = "LRT",
              parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_G_H1_H2<-results(desG,contrast = c("Population","HAL_I","HAL_II"),test = "Wald",
                    parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_G_H1_F<-results(desG,contrast = c("Population","HAL_I","FIL"),test = "Wald",
                   parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)


GL_G_H2_F<-results(desG,contrast = c("Population","HAL_II","FIL"),test = "Wald",
                   parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

#### E
x<-flist[[2]]
ddsE<- DESeqDataSet(se = se, design = as.formula(x[1]))
cat(paste(c("Start Running \"E Gene\" Analysis at",date(), "........\n")))
desE<-DESeq(ddsE,test="LRT", reduced= as.formula(x[2]), parallel = T)

resultsNames(desE)

GL_E<-results(desE,test = "LRT",
                  parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_E_H1_H2<-results(desE,contrast = c("SITE","BFL","KINGS"),test = "LRT",
              parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)


#### GXE
x<-flist[[3]]
ddsGE<- DESeqDataSet(se = se, design = as.formula(x[1]))
cat(paste(c("Start Running \"GE Gene\" Analysis at",date(), "........\n")))
desGE<-DESeq(ddsGE,test="LRT", reduced=  as.formula(x[2]), parallel = T)

resultsNames(desGE)

GL_GE<-results(desGE,test = "LRT",
                        parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_GE_B_K_H1_F<-results(desGE,name="PopulationHAL_I.SITEBFL",test = "Wald",
                        parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_GE_B_K_H2_F<-results(desGE,name="PopulationHAL_II.SITEBFL", test = "Wald",
                        parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

GL_GE_B_K_H1_H2<-results(desGE,contrast = list(c("PopulationHAL_I.SITEBFL",
                                                 "PopulationHAL_II.SITEBFL")),
                         test = "Wald",
                         parallel = T,alpha = 0.1,pAdjustMethod = "fdr",tidy = T)

save.image("Results/PH_Trans_Global.RData")
