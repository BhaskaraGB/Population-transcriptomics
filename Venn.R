
#install.packages("VennDiagramm")

library(VennDiagram)

#HALI_vs_HALII
GGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALI_vs_HALII.xlsx",sheetIndex = 1)
EGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALI_vs_HALII.xlsx",sheetIndex = 2)
GEGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALI_vs_HALII.xlsx",sheetIndex = 3)

venn.diagram(
x = list(GGene$row, EGene$row, GEGene$row),
  category.names = c("G" , "E" , "GXE"),
  filename = "HAL_I_vs_HAL_II.tiff",
  col = c("navyblue", "springgreen4","firebrick4"),
  fill = c("navyblue", "springgreen2","firebrick2"),
  alpha = 0.25,
  #main
  main= "HAL-I vs HAL-II",
  main.fontface="bold",
  main.fontfamily = "sans",
  main.cex = 1.5,
  main.pos = c(0.5,0.85),

  # Output features
  imagetype="tiff" ,
  height = 8 , 
  width = 8 , 
  resolution = 300,
  compression = "lzw",
  units = "in",
  
  #circle
  lwd = 2,
  lty = 1,
  
  #cat
  cex = 1.5,
  #fontfamily = "serif",
  fontface = "bold",
  cat.col = "black",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25, 180),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  rotation = 1,
  margin = 0.2
)


#HALI_vs_FIL
GGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALI_vs_FIL.xlsx",sheetIndex = 1)
EGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALI_vs_FIL.xlsx",sheetIndex = 2)
GEGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALI_vs_FIL.xlsx",sheetIndex = 3)

venn.diagram(
  x = list(GGene$row, EGene$row, GEGene$row),
  category.names = c("G" , "E" , "GXE"),
  filename = "HAL_I_vs_FIL.tiff",
  col = c("navyblue", "springgreen4","firebrick4"),
  fill = c("navyblue", "springgreen2","firebrick2"),
  alpha = 0.25,
  #main
  main= "HAL-I vs FIL",
  main.fontface="bold",
  main.fontfamily = "sans",
  main.cex = 1.5,
  main.pos = c(0.5,0.85),
  
  # Output features
  imagetype="tiff" ,
  height = 8 , 
  width = 8 , 
  resolution = 300,
  compression = "lzw",
  units = "in",
  
  #circle
  lwd = 2,
  lty = 1,
  
  #cat
  cex = 1.5,
  #fontfamily = "serif",
  fontface = "bold",
  cat.col = "black",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25, 180),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  rotation = 1,
  margin = 0.2
)

#HALII_vs_FIL
GGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALII_vs_FIL.xlsx",sheetIndex = 1)
EGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALII_vs_FIL.xlsx",sheetIndex = 2)
GEGene<-read.xlsx("PH_Trans_3Grp_DEG_cont_HALII_vs_FIL.xlsx",sheetIndex = 3)

venn.diagram(
  x = list(GGene$row, EGene$row, GEGene$row),
  category.names = c("G" , "E" , "GXE"),
  filename = "HAL_II_vs_FIL.tiff",
  col = c("navyblue", "springgreen4","firebrick4"),
  fill = c("navyblue", "springgreen2","firebrick2"),
  alpha = 0.25,
  #main
  main= "HAL-II vs FIL",
  main.fontface="bold",
  main.fontfamily = "sans",
  main.cex = 1.5,
  main.pos = c(0.5,0.85),
  
  # Output features
  imagetype="tiff" ,
  height = 8 , 
  width = 8 , 
  resolution = 300,
  compression = "lzw",
  units = "in",
  
  #circle
  lwd = 2,
  lty = 1,
  
  #cat
  cex = 1.5,
  #fontfamily = "serif",
  fontface = "bold",
  cat.col = "black",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25, 180),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  rotation = 1,
  margin = 0.2
)

