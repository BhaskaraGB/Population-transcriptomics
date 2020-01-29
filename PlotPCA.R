## Soon will release this as a package
## contact me if needed taslima@utexas.edu
## Wrutten by: Taslima Haque
## Jan 8th,2020

library(ade4)
library(factoextra)
library(magrittr)
library("FactoMineR")
library(tidyverse)


PlotPCA<-function (dat,
  axis1=1,
  axis2=2,
  shape=NULL,
  col.ind=NULL,
  ncp=5,
  mvar=NULL, # a vector of mapping factors which will be the same length for the number of variables that are used to construct PCA
  labmvar=NULL,
  colvar=NULL# a vector of color. Same lenth of mvar 
  ) {
  if (max(axis1,axis2) >5) {stop("axis number is greater than 5.\nSet ncp accordingly for your axis\n" )}
  if (length(mvar) != length(colvar)) {stop("length of mvar and colvar should be same")}
  if (!is.null(mvar) & is.null(labmvar)) {stop("mvar need a label for legend title")}
  
  axis<-c(axis1,axis2)
  seldat<-dat[,-c(shape,col.ind)]
  
  if ( !is.null(mvar) & dim(seldat)[2] != length(mvar)) {stop("Number of variables for PCA construction and the length of mvar not same. Can't map Variables")}
  
  res.pca <- PCA(seldat, graph = FALSE,ncp = ncp)
  perVar <- res.pca$eig[,2]
  ind <- data.frame(res.pca$ind$coord[, axis, drop = FALSE])
  var<-as.data.frame(res.pca$var$coord)
  origin=0
  origin <- rep(origin, nrow(var))
  dd <- cbind.data.frame(var, xstart = origin, ystart = origin)
  
  r <- min((max(ind[, 1]) - min(ind[, 1])/(max(var[, 1]) - min(var[, 1]))), (max(ind[, 2]) - min(ind[, 2])/(max(var[,2]) - min(var[, 2]))))
  scale = r * 0.7
  
  dd$xlab=dd$Dim.1*scale*1.08
  dd$ylab=dd$Dim.2*scale*1.08

  if ( !is.null(shape)) {ind$shape<-as.vector(dat[,shape])}
  if ( !is.null(col.ind)) {ind$col<-as.vector(dat[,col.ind])}
  xlab=paste("Axis",axis1," (",round(as.numeric(perVar[axis1]),2),"%)",sep = "")
  ylab=paste("Axis",axis2," (",round(as.numeric(perVar[axis2]),2),"%)",sep = "")
  
  floorOrceiling<-function(x) {
    if (x>= 0) {y<-ceiling(x); return(y)} else {
      y<-floor(x); return(y)
    }
  }
  
  MinMax<-function (ind) {
    lims<-c()
    for (i in 1:2) {
      if ( abs(min(ind[,i])) >= abs(max(ind[,i])) ) 
        { 
        vmax<-abs(floorOrceiling(min(ind[,i])) )
        vmin<-floorOrceiling(min(ind[,i])) 
        } else {
        vmax<-abs(floorOrceiling(max(ind[,i])) )
        vmin<-(floorOrceiling(max(ind[,i]))*-1)
        }
      lims<-append(lims,c(vmin,vmax))
    }
   return(lims) 
  }
  
  lims<-MinMax(ind)
  #shapelist<-c(8,1,6,17,0,2,3,4)
  ## Need to work on shape a bit
  shapelist<-c(16,17,18,19,20)
  ## And color too
  #collist<- c("#000000", "#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
  collist<- c( "#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
  ######### Plot individulas
  title<-paste("PCA on Axis ",axis1, " and ",axis2,sep = "")
  if (!is.null(shape) & !is.null(col.ind)){
    p<-ggplot(ind)+geom_point(aes(x=get(colnames(ind)[1]),y=get(colnames(ind)[2]),shape=ind$shape,colour=ind$col),size=3, alpha=0.75)+
      geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
      theme_bw() + labs(x =xlab, y=ylab,title = title,shape=as.name(colnames(dat)[shape]),colour=as.name(colnames(dat)[col.ind]))+ 
      scale_shape_manual(values=shapelist)+
      scale_color_manual(values=collist)+
      scale_x_continuous(limits = lims[1:2])+
      scale_y_continuous(limits = lims[3:4])+
      theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
            legend.title = element_text(size=10, face="bold"),
            legend.text = element_text(size=10),
            legend.position=c(1.05, 0.75),
            axis.title = element_text(size=12,face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            #panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text=element_text(size=12,vjust=1),
            plot.margin = margin(0.25, 5, 0.25, 0.25, "cm")
            #axis.ticks.x = element_blank()
      )
  }
  
  if (!is.null(shape) & is.null(col.ind)){
  p<-ggplot(ind)+geom_point(aes(x=get(colnames(ind)[1]),y=get(colnames(ind)[2]),shape=ind$shape),size=3, colour="black",alpha=0.75)+
    geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
    theme_bw() + labs(x =xlab, y=ylab,title = title,shape=as.name(colnames(dat)[shape]))+ 
    scale_shape_manual(values=shapelist)+
    scale_x_continuous(limits = lims[1:2])+
    scale_y_continuous(limits = lims[3:4])+
    theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10),
          legend.position=c(1.05, 0.75),
          axis.title = element_text(size=12,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12,vjust=1),
          plot.margin = margin(0.25, 5, 0.25, 0.25, "cm")
          #axis.ticks.x = element_blank()
    )
  }
  
  if (is.null(shape) & !is.null(col.ind)){
    p<-ggplot(ind)+geom_point(aes(x=get(colnames(ind)[1]),y=get(colnames(ind)[2]),colour=ind$col),size=3,alpha=0.75)+
      geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
      theme_bw() + labs(x =xlab, y=ylab,title = title,colour=as.name(colnames(dat)[col.ind]))+ 
      scale_color_manual(values=collist)+
      scale_x_continuous(limits = lims[1:2])+
      scale_y_continuous(limits = lims[3:4])+
      theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
            legend.title = element_text(size=10, face="bold"),
            legend.text = element_text(size=10),
            legend.position=c(1.05, 0.75),
            axis.title = element_text(size=12,face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            #panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text=element_text(size=12,vjust=1),
            plot.margin = margin(0.25, 5, 0.25, 0.25, "cm")
            #axis.ticks.x = element_blank()
      )
  }
  #### Plot Variables
  

  if(!is.null(mvar)) {
    p<-p+geom_segment(data=dd,aes(x=dd$xstart,y=dd$ystart,xend=dd$Dim.1*scale,yend=dd$Dim.2*scale),size=0.5,
                    arrow = arrow(length = unit(0.25, "cm")))+
    geom_text(data=dd, aes(x=dd$xlab, y=dd$ylab, label=dd$name), size=3,colour=col,fontface="bold")+
    labs(colour = "Treatment Stage")
    
  }
  return(p)
  }
