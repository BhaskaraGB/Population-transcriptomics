setwd("/home/taslima/data/JuengerLab/ALL_RNASeq/Bhashkar_RNASeq/V4")

rm(list=ls())

library(DESeq2)
library(ggplot2)
library(devtools)
library("BiocParallel")
library("vsn")
library(edgeR)
library(xlsx)

register(MulticoreParam(workers=8))

Ortho<-read.csv("tableS3_orthologDatabase.csv",check.names = F)
dim(Ortho)
Ortho<-na.omit(Ortho)
dim(Ortho)
Ortho$ID<-sapply(1:dim(Ortho)[1], function(x) paste( c("Ortho",sprintf("%06d", x)),collapse = "_"))

datH<-read.csv("Counts_HAL_rename.csv",check.names = F,row.names =1)
colnames(datH)<-gsub( "_-[0-9][0-9]|--[0-9][0-9]|_-[0-9]","",colnames(datH))
rownames(datH)<-gsub(".v2.1","",rownames(datH)[1:nrow(datH)])
datH<-datH[Ortho$HAL2,]

Design<-read.csv("Metasheet_lib.csv",row.names = 7,check.names = F,na.strings = "NA")
Design<-Design[-which(is.na(Design$location)),]
Design<-Design[-which(Design$SITE=="BFL_2"),]
Design<-Design[-which(rownames(Design) %in% c("PD19-30","PD19-5")),]

ID<-intersect(rownames(Design), colnames(datH))

dat<-datH[,ID]

Design<-Design[which(rownames(Design) %in% ID),]
dat<-dat[,rownames(Design)]
dat<-dat[-which(as.vector(rowMeans(dat,na.rm = F)) <1),]

dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = Design,
                              design = ~ 1  )

dds = estimateSizeFactors(dds )
sizeFactors(dds)

# Variance Stabilizing Transformation
vst<-vst(dds,blind = T)
dat<-as.data.frame(assay(vst))
tdat<-as.data.frame(t(dat))
head(tdat)
dim(tdat)
tdat$lib_ID<-rownames(tdat)

PheDes<-read.csv("GB_Metalib_Pheno.csv",check.names = F)
PheDes$lib_ID<-as.character(PheDes$lib_ID)
PheDes<-PheDes[-which(PheDes$location == "Mexico"),] #Remove Mexican samples
HAL2remove<-setdiff(unique(PheDes$genotype[which(grepl("HAL",PheDes$genotype))]),c("HAL2"))
FIL2remove<-setdiff(unique(PheDes$genotype[which(grepl("FIL",PheDes$genotype))]),c("FIL2"))
PheDes<-PheDes[-which(PheDes$genotype %in% c("GIL42A",HAL2remove,FIL2remove)),] # Remove outlier
PheDes$location[which(PheDes$location =="Coastal")] <-"Filipes"
PheDes$Population<-ifelse( PheDes$location %in% c("CenTex","West"),"HAL_I",ifelse(PheDes$location %in% c("Sym","TexMex"),"HAL_II","FIL"))
PheDes$Population<-factor(PheDes$Population,levels = c("FIL","HAL_I","HAL_II"))
PheExp<-merge(PheDes,tdat,by="lib_ID",all.x=F, all.y=F)
PheExp<-PheExp[,-c(1:2,4,7:15)]

ExpCount <- as.tibble(PheExp[,c(1:5)])  %>% group_by(genotype,location,SITE,Group,Population) %>% summarise(n=n())

GenCorPheExp <- as.tibble(PheExp[-which(ExpCount$n < 3),])  %>% group_by(genotype,location,SITE,Group,Population) %>% summarise_all(funs(mean))
rownames(GenCorPheExp)<-paste(GenCorPheExp$genotype,GenCorPheExp$SITE,sep = "_")
dim(GenCorPheExp)
tGenCorPheExp<-as.data.frame(t(GenCorPheExp[,-c(1:5)]))
colnames(tGenCorPheExp)<-rownames(GenCorPheExp)

#Get the variances of each gene to select top 1000 genes for PCA
var<-as.data.frame(apply(tGenCorPheExp,1,var,na.rm=T))
colnames(var)<-"GeneVariance"
#TopVarGene<-rownames(var)[order(var$GeneVariance,decreasing = T)][1:dim(var)[1]]
TopVarGene<-rownames(var)[order(var$GeneVariance,decreasing = T)][1:1000]
dat_sel<-tGenCorPheExp[which(rownames(tGenCorPheExp) %in% TopVarGene  ),]
colnames(dat_sel)<-paste(GenCorPheExp$genotype,GenCorPheExp$SITE,sep = "_")
dat<-cbind(GenCorPheExp[,c(5,3)],na.omit(t(dat_sel)))

p1<-PlotPCA(dat,axis1 = 1,axis2 = 2,shape = 2,col.ind=1,ncp = 5)
p2<-PlotPCA(dat,axis1 = 2,axis2 = 3,shape = 2,col.ind=1,ncp = 5)
p3<-PlotPCA(dat,axis1 = 3,axis2 = 4,shape = 2,col.ind=1,ncp = 5)
p4<-PlotPCA(dat,axis1 = 4,axis2 = 5,shape = 2,col.ind=1,ncp = 5)

grid.arrange(p1,p2,p3,p4, ncol=2, nrow = 2, top= "")

########### COV

PheDes<-read.csv("GB_Metalib_Pheno.csv",check.names = F)
PheDes$lib_ID<-as.character(PheDes$lib_ID)
PheDes<-PheDes[-which(PheDes$location == "Mexico"),] #Remove Mexican samples
PheDes<-PheDes[-which(PheDes$genotype %in% c("FIL216","GIL42A")),] # Remove outlier
PheDes$location[which(PheDes$location =="Coastal")] <-"Filipes"
PheDes$Population<-ifelse( PheDes$location %in% c("CenTex","West"),"HAL_I",ifelse(PheDes$location %in% c("Sym","TexMex"),"HAL_II","FIL"))
PheDes$Population<-factor(PheDes$Population,levels = c("FIL","HAL_I","HAL_II"))
PheDes$Group<-paste(PheDes$Population,PheDes$SITE,sep = "_")
PheExp<-merge(PheDes,tdat,by="lib_ID",all.x=F, all.y=F)
PheExp<-PheExp[,-c(1:2,4,6:15)]

ExpCount <- as.tibble(PheExp[,c(1:4)])  %>% group_by(genotype,SITE,Group) %>% summarise(n=n())

ExpCount<-ExpCount[-which(ExpCount$n < 3),] 
ExpCount$GenoSite<-paste(ExpCount$genotype,ExpCount$SITE,sep = "_")

PheDes$GenoSite<-paste(PheDes$genotype,PheDes$SITE,sep="_")
PheDes<-PheDes[which(PheDes$GenoSite %in% ExpCount$GenoSite ),]
rownames(PheDes)<-PheDes$lib_ID

CV<-function(x){ 
  cv<-round( 100*((sd(x,na.rm = T))/mean(x,na.rm=T)), 2)
  return(cv)
}
datCV<-PheExp %>% group_by(Group) %>% summarise_at(colnames(PheExp)[-(1:4)], CV)
df<-datCV %>% gather(GeneName, CV, colnames(datCV)[-1])
df$SITE<-rep("BFL",dim(df)[1])

datCVSite<-PheExp %>% group_by(SITE) %>% summarise_at(colnames(PheExp)[-(1:4)], CV)
dfSite<-datCVSite %>% gather(GeneName, CV, colnames(datCV)[-1])

ggplot() + 
  geom_line( data=df, mapping=aes(x=CV,color=Group,linetype=SITE),stat="density",size=0.5,alpha=1,show.legend  = T)+
  geom_line( data = dfSite,mapping = aes(x=CV,linetype=SITE),stat="density",size=0.8,alpha=0.75,show.legend = F)+
  theme_classic(base_size = 14)+labs(x="Coefficient of Variantion (CV)")+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                     #c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                     #c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                     #c("#000000", "#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7"),
                     name="GroupXSITE",
                     breaks=levels(factor(df$Group)),
                     labels= levels(factor(df$Group)))

