setwd("/home/taslima/data/JuengerLab/ALL_RNASeq/Bhashkar_RNASeq/V5")

rm(list=ls())

library(DESeq2)
library(ggplot2)
library(devtools)
library("BiocParallel")
library("vsn")
library(edgeR)
library(xlsx)
library(tidyverse)
library(FactoMineR)

register(MulticoreParam(workers=8))

Design<-read.csv("PHTrans_Final_Meta.csv",row.names = 6)
dat<-read.table("../V5/PHTrans_Count_HALREF.tab",header = T,stringsAsFactors = F,check.names = F)

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

PheDes<-Design[,-c(3)]
PheDes$lib_ID<-rownames(PheDes)
PheExp<-merge(PheDes,tdat,by="lib_ID",all.x=F, all.y=F)

ExpCount <- as.tibble(PheExp) %>% select(genotype,location,SITE,Group,Population)  %>% group_by(genotype,location,SITE,Group,Population) %>% summarise(n=n())

GenCorPheExp <- as.tibble(PheExp[-which(ExpCount$n < 3),])  %>% group_by(genotype,location,SITE,Group,Population) %>% summarise_all(funs(mean))
rownames(GenCorPheExp)<-paste(GenCorPheExp$genotype,GenCorPheExp$SITE,sep = "_")
dim(GenCorPheExp)
tGenCorPheExp<-as.data.frame(t(GenCorPheExp[,-c(1:19)]))
colnames(tGenCorPheExp)<-rownames(GenCorPheExp)

#Get the variances of each gene to select top 1000 genes for PCA
var<-as.data.frame(apply(tGenCorPheExp,1,var,na.rm=T))
colnames(var)<-"GeneVariance"
#TopVarGene<-rownames(var)[order(var$GeneVariance,decreasing = T)][1:dim(var)[1]]
TopVarGene<-rownames(var)[order(var$GeneVariance,decreasing = T)][1:1000]
dat_sel<-tGenCorPheExp[which(rownames(tGenCorPheExp) %in% TopVarGene  ),]
colnames(dat_sel)<-paste(GenCorPheExp$genotype,GenCorPheExp$SITE,sep = "_")
Expdat<-cbind(GenCorPheExp[,c(5,3)],na.omit(t(dat_sel)))

####### CODE FOR PCA
### It is raw right now
## Soon will update as a package

### Axis 1 and 2
axis1=1
axis2=2
shape=2
col.ind=1
axis<-c(axis1,axis2)

seldat<-Expdat[,-c(shape,col.ind)]
res.pca <- PCA(seldat, graph = FALSE,ncp = 5)
perVar <- res.pca$eig[,2]
ind <- data.frame(res.pca$ind$coord[, axis, drop = FALSE])
var<-as.data.frame(res.pca$var$coord)

ind$shape<-as.vector(Expdat[,shape])
ind$col<-as.vector(Expdat[,col.ind])
xlab=paste("Axis",axis1," (",round(as.numeric(perVar[axis1]),2),"%)",sep = "")
ylab=paste("Axis",axis2," (",round(as.numeric(perVar[axis2]),2),"%)",sep = "")

### Set the scale manually
fivenum(ind$Dim.1)
fivenum(ind$Dim.2)


## Need to work on shape a bit
shapelist<-c(21,24)

## And color too
collist<- c("#CB4335", "#566573", "#3949AB")

######### Plot individulas
title<-paste("PCA on Axis ",axis1, " and ",axis2,sep = "")

p1<-ggplot(ind)+geom_point(aes(x=ind[,1],y=ind[,2],shape=ind$shape,fill=ind$col),size=3,stroke=1.05, alpha=0.75,show.legend = F)+
    geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
    theme_bw() + labs(x =xlab, y=ylab,title = title,shape=as.name(colnames(Expdat)[shape]),fill=as.name(colnames(Expdat)[col.ind]))+ 
    scale_shape_manual(values=shapelist)+
    scale_fill_manual(values=collist)+
    scale_x_continuous(breaks = c(-30,-20,-10,0,10,20,30),limits = c(-34,34))+
    scale_y_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25))+
    theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10),
          legend.position=c(1.05, 0.75),
          axis.title = element_text(size=12,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12,vjust=1),
          plot.margin = margin(0.15, 0.15, 0.15, 0.15, "cm") # Change margin here to add/remove extra margin
          #axis.ticks.x = element_blank()
    )

############## Axis 2 and 3

axis1=2
axis2=3
shape=2
col.ind=1
axis<-c(axis1,axis2)

seldat<-Expdat[,-c(shape,col.ind)]
res.pca <- PCA(seldat, graph = FALSE,ncp = 5)
perVar <- res.pca$eig[,2]
ind2 <- data.frame(res.pca$ind$coord[, axis, drop = FALSE])
var2<-as.data.frame(res.pca$var$coord)

ind2$shape<-as.vector(Expdat[,shape])
ind2$col<-as.vector(Expdat[,col.ind])
xlab2=paste("Axis",axis1," (",round(as.numeric(perVar[axis1]),2),"%)",sep = "")
ylab2=paste("Axis",axis2," (",round(as.numeric(perVar[axis2]),2),"%)",sep = "")

### Set the scale manually
fivenum(ind2$Dim.2)
fivenum(ind2$Dim.3)


## Need to work on shape a bit
shapelist<-c(21,24)

## And color too
collist<- c("#CB4335", "#566573", "#3949AB")

######### Plot individulas
title<-paste("PCA on Axis ",axis1, " and ",axis2,sep = "")

p2<-ggplot(ind2)+geom_point(aes(x=ind2[,1],y=ind2[,2],shape=ind2$shape,fill=ind2$col),size=3,stroke=1.05, alpha=0.75,show.legend = F)+
  geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
  theme_bw() + labs(x =xlab2, y=ylab2,title = title,shape=as.name(colnames(Expdat)[shape]),fill=as.name(colnames(Expdat)[col.ind]))+ 
  scale_shape_manual(values=shapelist)+
  scale_fill_manual(values=collist)+
  scale_x_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25))+
  scale_y_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25))+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.position=c(1.05, 0.75),
        axis.title = element_text(size=12,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1),
        plot.margin = margin(0.15, 0.15, 0.15, 0.15, "cm") # Change margin here to add/remove extra margin
        #axis.ticks.x = element_blank()
  )

######## Get the Legend

## dummy plot

dum1<-ggplot(ind2)+geom_point(aes(x=ind2[,1],y=ind2[,2],shape=ind2$shape),size=3, stroke=1.25,alpha=1,show.legend = T)+
  geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
  theme_bw() + labs(x =xlab2, y=ylab2,title = title,shape=as.name(colnames(Expdat)[shape]))+ 
  scale_shape_manual(values=shapelist)+
  scale_x_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25))+
  scale_y_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25))+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        legend.position=c(1.05, 0.75),
        axis.title = element_text(size=12,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1),
        plot.margin = margin(0.25, 5, 0.25, 0.25, "cm")
        #axis.ticks.x = element_blank()
  )

dum2<-ggplot(ind2)+geom_point(aes(x=ind2[,1],y=ind2[,2],color=ind2$col),size=5,alpha=0.75,show.legend = T)+
  geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
  theme_bw() + labs(x =xlab2, y=ylab2,title = title,color=as.name(colnames(Expdat)[col.ind]))+ 
  scale_color_manual(values=collist)+
  scale_x_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25))+
  scale_y_continuous(breaks = c(-20,-10,0,10,20),limits = c(-25,25))+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        legend.position=c(1.05, 0.75),
        axis.title = element_text(size=12,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1),
        plot.margin = margin(0.25, 5, 0.25, 0.25, "cm")
        #axis.ticks.x = element_blank()
  )


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legendS <- get_legend(dum1)
legendC <- get_legend(dum2)

library(cowplot)

plot.with.inset <-
  ggdraw() +
  draw_plot(p1,x = 0, y = 0,width = 0.45,height = 1) +
  draw_plot(p2,x = 0.45, y = 0,width = 0.45,height = 1) +
  draw_plot(legendS, x = 0.62, y = .68, width = .3, height = .3)+
  draw_plot(legendC, x = 0.62, y = .55, width = .3, height = .3)

tiff("test.tiff",width=14,height=9,units="in",res=300)
plot.with.inset
dev.off()

pdf("test.pdf",width=14,height=9)
plot.with.inset
dev.off()

####### Size factor distribution

SF<-as.data.frame(sizeFactors(dds))
SF$lib_ID<-rownames(SF)
colnames(SF)[1]<-"sizeFactors"
PheSF<-PheDes[,c(2,13,15)]
SFdat<-merge(PheSF,SF,by="lib_ID")
ggplot(dat, aes(sizeFactors))+geom_density(aes(colour=Group,fill=Group),alpha=0.4)

ggplot() + 
  geom_line( data=SFdat, mapping=aes(x=sizeFactors,color=Group,linetype=SITE),stat="density",size=0.5,alpha=1,show.legend  = T)+
  geom_line( data = SFdat,mapping = aes(x=sizeFactors,linetype=SITE),stat="density",size=0.8,alpha=0.75,show.legend = F)+
  theme_classic(base_size = 14)+labs(x="Size Factors of libraries")+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                     name="GroupXSITE",
                     breaks=levels(factor(SFdat$Group)),
                     labels= levels(factor(SFdat$Group)))


########### COV
head(PheExp)[c(1:5),c(1:20)]
PheExp2<-PheExp[,-c(1,2,4:13)]

CV<-function(x){ 
  cv<-round( 100*((sd(x,na.rm = T))/mean(x,na.rm=T)), 2)
  return(cv)
}
datCV<-PheExp2 %>% group_by(Group) %>% summarise_at(colnames(PheExp2)[-(1:3)], CV)
df<-datCV %>% gather(GeneName, CV, colnames(datCV)[-1])
df$SITE<-rep("BFL",dim(df)[1])

datCVSite<-PheExp2 %>% group_by(SITE) %>% summarise_at(colnames(PheExp2)[-(1:3)], CV)
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


########## Pheno PCA

############ genetic correlation of pheno
head(PheDes)
PheDes<-PheDes[,c(1,2,14,4:12)]

#colnames(PheMat)<-c("PLOT_GL", "BIOMASS",  "TC",  "AREA","SLA", "FRESH", "DRY","RWC",  "OP", "MDWP")


### Change in to Numeric value
for (i in c(4:12)) {
  PheDes[,i]<-as.numeric(as.character(PheDes[,i]))
}

### Mask any RWC > 100 as NA
PheDes$RWC[which(PheDes$RWC >100)]<-NA

#write.csv(PheDes,"GB_Meta_647_withPhe.csv",row.names = F)

GenCorPhe <- as.tibble(PheDes)  %>% group_by(genotype,SITE,Population) %>% summarise_all(funs(mean))

PCAPhe<-GenCorPhe %>% ungroup() %>%  select(SITE,Population, BIOMASS,TC,AREA,SLA,FRESH,DRY,RWC,OP,MDWP) 
PCAPhe<-na.omit(PCAPhe)

### Distribution of pheno for two sites
PCAPhe %>%
  #keep(is.numeric) %>%                     # Keep only numeric columns
  #group_by(SITE) %>% 
  gather(key, value,BIOMASS,TC,AREA,SLA,FRESH,DRY,RWC,OP,MDWP) %>%                             # Convert to key-value pairs
  ggplot(aes(value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_density(aes(colour=SITE,fill=SITE),alpha=0.4) +
  #geom_line(aes(colour=SITE,fill=SITE),stat = "density",alpha=0.75) +
  scale_fill_manual(values = c("orange2","dodgerblue"))+
  scale_color_manual(values = c("orange3","navyblue"))+
  labs(x="Phenotype",y="Density")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12,face="bold"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1),
        strip.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(size=12,face = "bold")
        #axis.ticks.x = element_blank()
  )

### TAKE HOME MSG: AREA, DRY FREASH almost have the same distributions: may be they all are leaf traits therefore have the same pattern
### MDWP and OP totally bimodal if we consider total population; this bimodal pattern is for SITE
### sd of RWC of KINGS are higher that BFL therefore have very different disribution and not has a clear bimodal separartion


####### RUN PCA on genetically identical pooled samples
####### Will change this PCA later
###### SITE should be in different color
library(ade4)
library(factoextra)
library(magrittr)
library(gridExtra)


res.pca <- dudi.pca(PCAPhe[,-c(1:2)],
                    scannf = FALSE,   # Hide scree plot
                    nf = 5 ,center = T,scale = T           # Number of components kept in the results
)

PCAPhe$Group<-paste(PCAPhe$Population,PCAPhe$SITE,sep = "_")
fviz_pca_biplot(res.pca, label ="var",axes = c(1,2),
                geom.var = c("arrow","text"),
                geom.ind = "point",
                col.var = "orangered2",
                pointshape = 21, 
                pointsize = 2.5,
                invisible = "quali",
                palette= c( "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"),
                #palette = colorRampPalette(brewer.pal(name="Paired", n = 16))(16),
                #palette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                #            "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
                fill.ind=PCAPhe$Group,
                legend.title = list(fill = "Population:Site"),
                #select.var = list(contrib = 10),
                repel = TRUE,
                addEllipses=F, ellipse.level=0.75)+labs(title="PCA on Axis 1 and 2")+
  theme(plot.title = element_text(hjust = 0.5))  

