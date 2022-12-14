library(topGO)

RunTopGO<-function(GOfile,info,geneList) {

geneID2GO <- readMappings(file = GOfile)
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
allRes_BP<-allRes_BP[which(allRes_BP$padj < as.vector(quantile(allRes_BP$padj,prob=0.1))),] # 10% FDR
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

return(rbind(allRes_BP,allRes_CC,allRes_MF))
}
