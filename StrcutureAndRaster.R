setwd("/home/taslima/data/JuengerLab/ALL_RNASeq/Bhashkar_RNASeq/V5")
########## PH Plot
library(tidyverse)
library(ggplot2)
##### Structure

Struc<-read.csv("Structure_for_Raster.csv")

Struc<-Struc[-which(is.na(Struc$Population)),]

#swap K2 and K3 from structure output

Struc<-as.tibble(Struc) %>% arrange(desc(K1,K3,K2)) %>% mutate( Pop= case_when(
  K1 >= 0.55 ~ "K1",
  K2 >=0.55 ~ "K2", 
  K3 >= 0.55 ~ "K3", 
  TRUE ~ "Mixed"
  
)) %>% arrange (Pop) %>% mutate(Index = seq(1:dim(Struc)[1]))



(str<-as.tibble(Struc) %>% gather(K1,K2,K3,key = "K",value = "Prob") %>% arrange(Index) %>% 
    ggplot(aes(x=Index,y=Prob,fill=K))+geom_bar(width = 1,position="stack", stat="identity",show.legend = F)+
    scale_fill_manual(values = c("#566573", "#3949AB","#CB4335"),labels=c("HAL-I","HAL-II","FIL"))+
    labs(x="",y="",fill="Population Group")+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA))
)

####### Plot in Map
library(raster)
library(maps)
library(rgdal)
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
library(scatterpie)


fivenum(Struc$Latitude)
fivenum(Struc$Longitude)

ext<-extent(-115, -95, 24, 38)

for (i in 1:12) {
  name<-paste("prec",i,sep="_")
  file=paste("/home/taslima/data/JuengerLab/RTeaching/OurCodingClub/SpatialR/wc2.0_2.5m_prec/wc2.0_2.5m_prec_",sprintf("%02d",i),".tif",sep="")
  cat(paste("reading \"",file ,"\" ......\n",sep = ""))
  tmp<-raster(file)
  tmp<-crop(tmp,ext)
  assign(name, tmp)
}

prec<-mean(prec_1,prec_2,prec_3,prec_4,prec_5,prec_6,prec_7,prec_8,prec_9,prec_10,prec_11,prec_12)


USA <- getData('GADM', country="USA", level=1)
prec<-crop(prec,bbox(USA))
prec<-mask(prec,USA)
prec2<-as.data.frame(prec,xy=T)
crs(prec)
prec3<-prec2[-which(is.na(prec2$layer)),]

mapdata <- read.csv("state-medal-count.csv", header=TRUE, stringsAsFactors=FALSE)
states <- map_data("state")
substates<-states[which(states$region %in% c("texas","new mexico","arizona","oklahoma")),]


p<-ggplot() +
  geom_point(data = prec3 , aes(x = x, y = y, color = layer),alpha=0.65) +
  scale_color_gradientn(colours =  terrain.colors(12))+
  geom_polygon(data=substates,aes(x = long, y = lat,group = group),alpha=0, color = "white") + 
  coord_quickmap()+
  labs(x="",y="",color="Precipitation",fill="SubPopulation")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))


Struc$K1<-Struc$K1*100
Struc$K2<-Struc$K2*100
Struc$K3<-Struc$K3*100

Struc$radius <-rep(0.25,dim(Struc)[1])

q<-p+ geom_scatterpie(aes(x=Longitude, y=Latitude,r=radius),data = Struc, cols=c("K1","K2","K3"), color=NA,alpha=0.75)+
  scale_fill_manual(values = c("#566573", "#3949AB","#CB4335"),labels=c("HAL-I","HAL-II","FIL"))+
  coord_quickmap()

mapdata <- read.csv("state-medal-count.csv", header=TRUE, stringsAsFactors=FALSE)
states <- map_data("state")
modstates<-states
modstates$Flag<-rep(NA,dim(modstates)[1])
modstates$Flag[which(modstates$region %in% c("texas","new mexico","arizona","oklahoma"))]<-"Sel"

sm<-ggplot(modstates)+geom_polygon(aes(x=long,y=lat,group=group,fill=Flag),alpha=0.5,color="grey30",size=0.2,show.legend = F)+
  scale_fill_manual(values = c("grey40","white"))+#theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA))+
  labs(x="",y="")+coord_quickmap()


library(cowplot)

plot.with.inset <-
  ggdraw() +
  draw_plot(q,x = 0, y = 0) +
  draw_plot(sm, x = 0.05, y = 0.2, width = .25, height = .25)+
  draw_plot(str, x = 0.01, y = 0.01, width = .6, height = .2)

tiff("PopStrcut.tiff",width=14,height=9,units="in",res=300)
plot.with.inset
dev.off()

