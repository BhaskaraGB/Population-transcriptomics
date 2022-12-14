setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/PHtrans/")

rm(list=ls())

library(qpctools)
library(viridis)
library(qvalue)
library(tidyverse)
library(ggplot2)

############ Get Count Data and make mean center count

Design<-read.csv("Data/PHTrans_Final_Meta_548lib_FullStat.csv",row.names = 1)
dat<-read.table("Data/PHTrans_Count_HALREF_ALL.tab",header = T,stringsAsFactors = F,check.names = F)

dat<-dat[,rownames(Design)]
dat<-dat[-which(as.vector(rowMeans(dat,na.rm = F)) <1),]

#Detected genes
#DetectedGene<-read.csv("Data/PHTrans_Expressed_GeneID.csv")

## Filter genes present only high REsites
RE<-read.delim("DBs/PhalliiHAL_496_v2.1.gene.REsites.bed",header = F)
colnames(RE)<-c("Chr","Start","End", "Score","Strand","GeneID")

dat<-dat[which(rownames(dat) %in% RE$GeneID ),]

#### Now rename genotype as 3-letter code
# Design$genotype<-substr(Design$genotype,1,3)
Design<-Design[-which(Design$genotype %in% c("COR","COL","OPV", "BUR","KEN","LME","SLS","SMF")),]
Design<-Design[,c('SITE','genotype','Population')]
Design$LIBID<-rownames(Design)

#### Merge with Expression count
library(DESeq2)
dat<-dat[,rownames(Design)]

dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = Design,
                              design = ~ 1  )

dds = estimateSizeFactors(dds )
SizeFac<-as.vector(sizeFactors(dds))


tdatSFCen<-as.data.frame(t(datSF))
tdatSFCen$LIBID<-rownames(tdatSFCen)

###### OR USE VST
vst<-vst(dds,blind = T)
datSF<-as.data.frame(assay(vst))
tdatSFCen<-as.data.frame(t(datSF))
tdatSFCen$LIBID<-rownames(tdatSFCen)

require(data.table)
SmallDes<-Design[,c('SITE','genotype','LIBID')]
datExp<-data.table(merge(SmallDes,tdatSFCen,by="LIBID"))


### helpful: https://riptutorial.com/data-table/example/13084/using--sd-and--sdcols

## Average of genotypes irrespective of SITE
datExpAvg<-datExp[order(genotype), lapply(.SD,mean), by = genotype, .SDcols = -c('LIBID','SITE')]

## Average of genotypeXSITE
datExpAvgSite<-datExp[order(genotype,SITE), lapply(.SD,mean), by = .(genotype,SITE), .SDcols = -c('LIBID')]

## Average of Average (this is better, unbiased for the sample size on each site)
datExpAvgAvg<-datExpAvgSite[order(genotype,SITE), lapply(.SD,mean), by = genotype,.SDcols = -c('SITE')]

## Difference between SITES: KINGS-BFL
datExpDiffSite<-datExpAvgSite[order(genotype,SITE), lapply(.SD,diff), by = genotype,.SDcols = -c('SITE')]

##get the needed Meta information
Meta<-as.data.frame(unique(Design[,c('genotype','Population')]))
rownames(Meta)<-Meta$genotype


##### Build Kinship
VCF<-read.delim("Data/PHNatAcc_AllChr_SNP_m90_WithHAL_PH86_DP3_mac3_ExcHET_50K_RE_mat.vcf",header = T,sep = "\t",check.names = T)

colnames(VCF)
colnames(VCF)[1]<-"CHROM"

colnames(VCF)[10:ncol(VCF)]<-substr(colnames(VCF)[10:ncol(VCF)],1,3)
#rownames(VCF)<-paste(VCF$CHROM,VCF$POS,sep = "_")
rownames(VCF)<-seq(1:nrow(VCF))


#Gen file
myG<-t(VCF[,-c(1:9)])

myG<-myG[ , apply(myG, 2, function(x) !any(is.na(x)))]
#myG<-myG[,-which(colMeans(myG) ==1)]
myG<-as.matrix(myG)

make_k_complete <- function(myG){
  scaleFactor = 1/sqrt(colMeans(myG) * (1 - colMeans(myG)))
  myS = matrix(0, nrow = dim(myG)[2], ncol = dim(myG)[2])
  diag(myS) = scaleFactor
  myM = dim(myG)[1]
  myT = matrix(data = -1/myM, nrow = myM - 1, ncol = myM)
  diag(myT) = (myM - 1)/myM
  myGstand = myT %*% myG %*% myS
  myK = cov(t(myGstand))
  return(myK)
}

K<-make_k_complete(myG)
myK<-K
row.names(myK) = rownames(myG)[-c(dim(myG)[1])]
colnames(myK) = rownames(myG)[-c(dim(myG)[1])]

eigF = eigen(myK)
myU = eigF$vectors
myLambdas = eigF$values


dim(myK)
length(colnames(myK))

#calculate the PC cutoffs we're using for tests of selection
varexp = myLambdas/sum(myLambdas)
sumexp = sapply(1:length(varexp), function(x){sum(varexp[1:x])})

#get cutoffs for how many pcs to look at
pcmax = which(sumexp > 0.25)[1]
plot(myLambdas, bty="n", xlab = "PCs", ylab = "Eigenvalues")
abline(v = pcmax, col = viridis(6)[3], lwd=2)
abline(v = which(sumexp > 0.25)[1], col = viridis(6)[2], lwd=2)
abline(v = which(sumexp > 0.20)[1], col = viridis(6)[1], lwd=2)
abline(v = which(sumexp > 0.50)[1], col = viridis(6)[3], lwd=2)
abline(v = which(sumexp > 0.80)[1], col = viridis(6)[3], lwd=2)
abline(v = which(sumexp > 0.90)[1], col = viridis(6)[4], lwd=2)

getqvalues <- function(ptable){
  qobj = qvalue(p = c(ptable))
  myqvals = matrix(qobj$qvalues, nrow=dim(ptable)[1])
  return(myqvals)
}

calcQpc<-function (myZ, myU, myLambdas, myPCcutoff=myPCcutoff, tailCutoff = 0.9, vapcs =vapcs) 
{
  myTailCutoff = round(tailCutoff * length(myLambdas))
  pcm = which(sapply(1:length(myLambdas), function(x) {
    sum(myLambdas[1:x])/sum(myLambdas)
  }) > myPCcutoff)[1]
  myZ = myZ[1:dim(myU)[1]] - mean(myZ)
  myCm = (myZ %*% myU)/sqrt(myLambdas)
  myQm = sapply(1:pcm, function(n) {
    #var0(myCm[n])/var0(myCm[(vapcs+1):myTailCutoff])
    var0(myCm[n])/var0(myCm[(myTailCutoff - vapcs):myTailCutoff])
  })
  myPs = sapply(1:pcm, function(x) {
    pf(myQm[x], 1, vapcs, lower.tail = F)
  })
  retdf = list(cm = myCm, qm = myQm, pvals = myPs)
  return(retdf)
}

# pcm = which(sapply(1:length(myLambdas), function(x) {
#   sum(myLambdas[1:x])/sum(myLambdas)
# }) > 0.3)[1]

myPCcutoff=0.25 # 6 PCs
vapcs = 45  ## Almost 80% from the upper tail

pcm = which(sapply(1:length(myLambdas), function(x) {
  sum(myLambdas[1:x])/sum(myLambdas)
}) > myPCcutoff)[1]

############## REORDER META Meta[match(rownames(myK),Meta$genotype),]
Meta<-Meta[match(rownames(myK),Meta$genotype),]


##################### AVERAGE
#bluptableAvg<-as.data.frame(datExpAvg[match(rownames(myK),datExpAvg$genotype),])
bluptableAvg<-as.data.frame(datExpAvgAvg[match(rownames(myK),datExpAvgAvg$genotype),])
bluptableAvg<-merge(Meta,bluptableAvg,by="genotype",all.x=F,all.y=T)
rownames(bluptableAvg)<-bluptableAvg$genotype
bluptableAvg<-as.data.frame(bluptableAvg[match(rownames(myK),bluptableAvg$genotype),])

mydfsAvg = apply(bluptableAvg[,-c(1:3)], 2, function(x){calcQpc(
  myZ = x, myU = eigF$vectors, myLambdas = eigF$values, myPCcutoff=myPCcutoff, vapcs = vapcs, tailCutoff = 0.9
)})


allpvalsAvg = sapply(1:length(mydfsAvg), function(x){mydfsAvg[[x]]$pvals})
myqvalsAvg = getqvalues(allpvalsAvg)
length(which(myqvalsAvg < 0.1))

################### DIFFERENCE
bluptableDiff<-as.data.frame(datExpDiffSite[match(rownames(myK),datExpDiffSite$genotype),])
bluptableDiff<-merge(Meta,bluptableDiff,by="genotype",all.x=F,all.y=T)
rownames(bluptableDiff)<-bluptableDiff$genotype
bluptableDiff<-as.data.frame(bluptableDiff[match(rownames(myK),bluptableDiff$genotype),])

mydfsDiff = apply(bluptableDiff[,-c(1:3)], 2, function(x){calcQpc(
  myZ = x, myU = eigF$vectors, myLambdas = eigF$values, myPCcutoff=myPCcutoff, vapcs = vapcs, tailCutoff = 0.9
)})


allpvalsDiff = sapply(1:length(mydfsDiff), function(x){mydfsDiff[[x]]$pvals})
myqvalsDiff = getqvalues(allpvalsDiff)
length(which(myqvalsDiff < 0.1))



############ GATHER RESULTS and REFORMAT 

calcCIs <- function(myName, myBlups=bluptable, myU=eigF$vectors, myLambdas=eigF$values,vapcs=vapcs,tailCutoff=tailCutoff){
  myTailCutoff = round(tailCutoff * length(myLambdas))
  myZ = myBlups[,myName]
  myZ = myZ - mean(myZ)
  myBm = myZ %*% myU
  myCm = myBm/sqrt(myLambdas)
  #myVa = var0(myCm[(vapcs+1):length(myLambdas)])
  myVa = var0(myCm[(vapcs+1):myTailCutoff])
  myCI = sqrt(myVa*myLambdas)
  return(myCI)}


gather_QPC_result<-function(bluptable,pcm,eigV=eigF$values,pvals,qvals,FDRcut=FDRcut) {
  #print(FDRcut)
  #print(pcm)
  GeneIndex<-((which(unlist(qvals) < FDRcut)%/%pcm)+1+3) #+3 column of geno,pop,subpop
  PCIndex<-(which(unlist(qvals) < FDRcut)%%pcm)
  PCIndex[which(PCIndex==0)]<-pcm
  pvalsel<-unlist(pvals)[which(unlist(qvals) < FDRcut)]
  qvalsel<-unlist(qvals)[which(unlist(qvals) < FDRcut)]
  GeneList<-colnames(bluptable)[GeneIndex]
  GeneMatrix<-as.data.frame(cbind(GeneName=GeneList,PCIndex=PCIndex,Pvalue=pvalsel,Qvalue=qvalsel))
  return(GeneMatrix)
}

gather_QPC_result_pval<-function(bluptable,pcm,eigV=eigF$values,pvals,qvals,Pvalcut=Pvalcut) {
  #print(FDRcut)
  #print(pcm)
  GeneIndex<-((which(unlist(pvals) < Pvalcut)%/%pcm)+1+3) #+3 column of geno,pop,subpop
  PCIndex<-(which(unlist(pvals) < Pvalcut)%%pcm)
  PCIndex[which(PCIndex==0)]<-pcm
  pvalsel<-unlist(pvals)[which(unlist(pvals) < Pvalcut)]
  qvalsel<-unlist(qvals)[which(unlist(pvals) < Pvalcut)]
  GeneList<-colnames(bluptable)[GeneIndex]
  GeneMatrix<-as.data.frame(cbind(GeneName=GeneList,PCIndex=PCIndex,Pvalue=pvalsel,Qvalue=qvalsel))
  return(GeneMatrix)
}

GenMatdiff<-gather_QPC_result(bluptableDiff,pcm,eigV = eigF$values,allpvalsDiff,myqvalsDiff,FDRcut = 0.1)
GenMatdifflen<-gather_QPC_result_pval(bluptableDiff,pcm,eigV = eigF$values,allpvalsDiff,myqvalsDiff,Pvalcut = 0.05)


info<-read.delim("DBs/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")
info<-info[,c(2,11:13,16)]
info<-info[!duplicated(info$locusName),]
dim(info)

GENEDIFFtabAnnot<-merge(GenMatdiff,info,by.x="GeneName",by.y="locusName",all.x=T,all.y=F)
GENEDIFFtabAnnotlen<-merge(GenMatdifflen,info,by.x="GeneName",by.y="locusName",all.x=T,all.y=F)

write.csv(GENEDIFFtabAnnotlen,"Results/QPC_Expression_Divergent_DIFFLEN_RE_v9.csv",
          row.names = F)


GenMatAve<-gather_QPC_result(bluptableAvg,pcm,eigV = eigF$values,pvals = allpvalsAvg,qvals = myqvalsAvg,FDRcut = 0.1)
GENEAVGtabAnnot<-merge(GenMatAve,info,by.x="GeneName",by.y="locusName",all.x=T,all.y=F)


write.csv(GENEAVGtabAnnot,"Results/QPC_Expression_Divergent_AVG_v9_RE.csv",
          row.names = F)



######### Plot Specific gene
GenMat<-GenMatAve
bluptable<-bluptableAvg
PCMat<-as.data.frame(eigF$vectors[,c(1:pcm)])
rownames(PCMat)<-rownames(myK)
colnames(PCMat)<-paste("PC",seq(1:pcm),sep = "")
PCMat$genotype<-rownames(PCMat)
PCBlupMat<-merge(PCMat,bluptable,by="genotype")

i=c(325)

myName<-as.character(GenMat$GeneName[i[1]])
myPC=paste("PC",as.numeric(as.character(GenMat$PCIndex[i[1]])),sep = "")
PlotMat<-as.data.frame(PCBlupMat[,c('genotype','Population','subpopulation',myPC,myName)])
PlotMat$subpopulation<-as.character(PlotMat$subpopulation)
PlotMat$subpopulation<-factor(PlotMat$subpopulation)
colnames(PlotMat)[c(4,5)]<-c("PC","Expression")
Intercept=mean(PlotMat$Expression)
myCI<-calcCIs(myName, myBlups=PCBlupMat, myU=eigF$vectors, myLambdas=eigF$values,vapcs=vapcs,tailCutoff=0.9)[1]
Slope=1.96*myCI
cat("Gene Name=", myName,"Mean value= ",Intercept,"CI=",myCI,"\n")

p1<-ggplot(PlotMat,aes(x=PC,y=Expression,colour=Population))+geom_point(size=3,alpha=0.75,show.legend = F)+
  geom_smooth(method = "lm",colour="black",se = F,size=0.5,linetype=1,fill="transparent")+
  geom_abline(intercept = Intercept,slope = Slope,linetype=2,color=viridis(6)[1])+
  geom_abline(intercept = Intercept,slope = -Slope,linetype=2,color=viridis(6)[1])+
  scale_colour_manual(values = c(HAL_I = "#4040a1",HAL_II = "#bd5734", FIL = "#5b9aa0"))+
  labs(title=myName,x=myPC,y="Normalized expression")+
  theme_classic(base_size = 14)+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5,size=14, face="bold"),
        axis.text = element_text(size=12, face="bold",color="black"),
        axis.title = element_text(face="bold"))


########## RUN IND PC GO
library(readxl)
library(xlsx)

info<-read.delim("DBs/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")
info<-info[,c(2,13,16)]
info<-info[!duplicated(info$locusName),]
dim(info)

GODB<-"DBs/PH_HAL_v2_GO_subset.tab"


RunGO<-function(GODB=GODB,info=info,geneList=geneList,FDR=FDR){
  library(topGO)
  #library(Rgraphviz)
  
  #read customized annotation file
  geneID2GO <- readMappings(file = GODB)
  
  ## Test with random p-value
  geneNames <- info$locusName
  
  myInterestingGenes<-as.data.frame(geneList)
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
  #allRes_BP<-allRes_BP[which(allRes_BP$padj < as.vector(quantile(allRes_BP$padj,prob=0.1))),] # 1% FDR
  allRes_BP$Attribute<-rep("BP",dim(allRes_BP)[1])
  
  ## MF
  
  GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=10)
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  pvalFis <- score(resultFisher)
  geneData(resultFisher)
  allRes_MF <- GenTable(GOdata, classic = resultFisher, ranksOf = "classic",topNodes=as.vector(attr(resultFisher,"geneData"))[4],numChar=100)
  allRes_MF$padj<-p.adjust(as.numeric(allRes_MF$classic),method = "fdr")
  #allRes_MF<-allRes_MF[which(allRes_MF$padj < as.vector(quantile(allRes_MF$padj,prob=0.1))),]
  allRes_MF$Attribute<-rep("MF",dim(allRes_MF)[1])
  
  ##CC
  GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=10)
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.stat)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes_CC <- GenTable(GOdata, classic = resultFisher, ranksOf = "classic",topNodes=as.vector(attr(resultFisher,"geneData"))[4],numChar=100)
  allRes_CC$padj<-p.adjust(as.numeric(allRes_CC$classic),method = "fdr")
  #allRes_CC<-allRes_CC[which(allRes_MF$padj < as.vector(quantile(allRes_CC$padj,prob=0.1))),]
  allRes_CC$Attribute<-rep("CC",dim(allRes_CC)[1])
  
  MergeallRes<-rbind(allRes_BP,allRes_CC,allRes_MF)
  #MergeallRes<-MergeallRes[which(MergeallRes$padj <= FDR),]  
  return(MergeallRes)
  
}

GENEAVGtabAnnot<-read.csv("Results/QPC_Expression_Divergent_AVG_v9_RE.csv",stringsAsFactors = F)
PClev<-sort(as.numeric(unique(GENEAVGtabAnnot$PCIndex)))

for (i in PClev) {
  PCGen<-GENEAVGtabAnnot$GeneName[which(GENEAVGtabAnnot$PCIndex==i)]
  PCGO<-RunGO(GODB=GODB,info=info,geneList=PCGen,FDR=0.5)
  sheetName=paste("PC",i,sep = "_")
  if (i==1) {
    write.xlsx(PCGO,file = "Results/QPC_v9_RE_GO.xlsx",sheetName = sheetName,row.names = F,append = F) 
  } else {
    write.xlsx(PCGO,file = "Results/QPC_v9_RE_GO.xlsx",sheetName = sheetName,row.names = F,append = T) 
  }
  
}



