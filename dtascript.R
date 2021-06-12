#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.1.st.v1)
library(hugene11sttranscriptcluster.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr)
library(tidyr)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
#library(devtools)


sample.names = c(  "Adult10","Adult11","Adult2","Adult4","Adult6","Pedo1", 
                  "Pedo3", "Pedo5", "Pedo7", "Pedo8", "Pedo9"  )

list = list.files('C:/Users/tapan/JUpyter/Biopython/GSE58315_RAW/',full.names = TRUE)
raw_data <- oligo::read.celfiles(filenames = list,sampleNames = sample.names)
affy.raw <- raw_data


l = list("Adult","Adult","Adult","Adult","Adult","Pre","Pre","Pre","Pre","Pre","Pre")



arrayQualityMetrics(expressionset = affy.raw, outdir =
                      "QC_report_raw", force = TRUE, do.logtransform = TRUE)
########      PCA     ######################################
exp_raw <- log2(Biobase::exprs(raw_data))


PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
PCA_raw2 <- princomp(exp_raw)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

for (i in 1:29){
  dataGG <- data.frame(PC1 = PCA_raw$x[,i], PC2 = PCA_raw$x[,i+1],Phenotype = unlist(l))
  
  # dataGG$Disease <- l 
  #dataGG
  
  
  print(ggplot(dataGG, aes(PC1, PC2)) +
          geom_point(aes(colour = Phenotype))+
          ggtitle("PCA plot of the log-transformed raw expression data") +
          xlab(paste0(cat("PC-", i, ", VarExp: ",percentVar[i], "%"))) +
          ylab(paste0(cat("PC-", i+1, ", VarExp: ", percentVar[i+1], "%"))))
  
}


dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],Phenotype = unlist(l))
png(file = "PCA0.png", width = 3000, height
    = 4000, res = 600)
(ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = Phenotype, size=3))+
    ggtitle("PCA plot of the log-transformed raw expression data") +
    xlab(paste0("PC-1, VarExp: ",percentVar[1], "%")) +
    ylab(paste0("PC-2, VarExp: ", percentVar[2], "%")))

dev.off()

#########################################################################

###############   RMA   #########################

eset <- rma(affy.raw)
write.exprs(eset,file="rma_norm_expr.txt")

#######  PCA after RMA

exp <- log2(Biobase::exprs(eset))

PCA <- prcomp(t(exp), scale. = FALSE)
PCA_2 <- princomp(exp)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],Phenotype = unlist(l))
png(file = "PCA1.png", width = 3000, height
    = 4000, res = 600)
(ggplot(dataGG, aes(PC1, PC2)) +
    geom_point(aes(colour = Phenotype, size=3))+
    ggtitle("PCA plot of the normalized expression data") +
    xlab(paste0("PC-1, VarExp: ",percentVar[1], "%")) +
    ylab(paste0("PC-2, VarExp: ", percentVar[2], "%")))

dev.off()

############ Gene Annotation   #############################

keytypes(hugene11sttranscriptcluster.db)
gns <- AnnotationDbi::select(hugene11sttranscriptcluster.db, keys(hugene11sttranscriptcluster.db), c("ENTREZID", "SYMBOL","PROBEID"))

gns <- subset(gns, !is.na(SYMBOL))

gns
gns <- gns[!duplicated(gns[,1]),]
head(gns)
row.names(gns) = gns$PROBEID
gns = gns[,-1]
head(gns)


expr <- data.frame(exprs(eset))
expr
expr.anno <- merge(gns, expr, by.x=0, by.y=0, all=TRUE)
head(expr.anno)

################## QUALITY METRICS   #######################

arrayQualityMetrics(expressionset = eset, outdir =
                      "QC_report", force = TRUE, do.logtransform = TRUE)

################


##############  GENE FILTERING ###############################

library(genefilter)

eset.filt = varFilter(eset, var.func=IQR, var.cutoff
                      = 0.25, filterByQuantile = TRUE)
nrow(eset.filt)
head(eset.filt)


################# Design Matrix
sample.groups <- factor(c("Adult","Adult","Adult","Adult","Adult","Pre","Pre","Pre","Pre","Pre","Pre"),
                         levels = c("Adult", "Pre"))
design.mat <- model.matrix(~0 + ~sample.groups)
colnames(design.mat)<- c("Adult", "Pre")


#################  Create Contrast Matrix

contrast.matrix <- makeContrasts(Adult-Pre, Pre-Adult, levels=design.mat)

fit <- lmFit(eset.filt, design.mat)
fit<- contrasts.fit(fit, contrasts = contrast.matrix)

fit <- eBayes(fit)
head(fit)

fit$genes <- gns [row.names(gns) %in% row.names(fit$t),]
fit.table=topTable(fit, coef=1, adjust.method = "BH", number=nrow(fit))
write.table(fit.table, file="fit_1.csv", sep=',', quote=FALSE, row.names=FALSE, col.names=TRUE)

#topTable(fit, coef=1, adjust.method = "BH", sort.by="logFC",lfc = 1, p.value = 0.005, number=10)
#topTable(fit$coefficients,number = 10)


##### LIST OF DEGS
de.gene = topTable(fit, coef=1, adjust.method = "BH", sort.by="logFC", 
                   lfc = 1,  number=nrow(eset.filt))
de.gene <- subset(de.gene, !is.na(SYMBOL))

write.table(de.gene, file="de_gene.csv", sep=',' , quote=FALSE, row.names=FALSE, col.names=TRUE)
#### UP AND DOWN REGULATED GENES ############
head(de.gene)

up.gene = de.gene[which(de.gene$logFC > 0), ]
up.gene = up.gene[order(round(as.numeric(up.gene$logFC)),decreasing = TRUE, na.last = TRUE),]
head(up.gene[order(round(as.numeric(up.gene$logFC)),decreasing = TRUE, na.last = TRUE),])


down.gene = de.gene[which(de.gene$logFC < 0), ]



head(down.gene[order(round(as.numeric(down.gene$logFC)),decreasing = FALSE, na.last = TRUE),])
down.gene = down.gene[order(round(as.numeric(down.gene$logFC)),decreasing = FALSE, na.last = TRUE),]
write.table(up.gene, file="up_gene.csv", sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(down.gene, file="down_gene.csv", sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)





######## VOLCANO

png(file = "volcano_plot1.png", width = 4000, height = 3000, res = 600)
volcanoplot(fit, coef = 2, highlight = 20)
#points(de.gene$logFC, de.gene$B, cex=.5, col=2, pch=19)
dev.off()



library(pheatmap)
de.gene.expr = merge(de.gene, expr[row.names(expr) %in%
                                     row.names(de.gene),], by.x=0, by.y=0, all=TRUE)
png(file = "Fig.de_gene.heatmap.png", width = 3000, height
    = 4000, res = 600)
pheatmap(as.matrix(de.gene.expr[,10:17]), scale="row",
         labels_row = de.gene.expr$SYMBOL)
dev.off()


############################ GO 
up.gene = de.gene.expr[which(de.gene.expr$logFC > 0), ]
up.gene = up.gene[order(round(as.numeric(up.gene$logFC)),decreasing = TRUE, na.last = TRUE),]
head(up.gene[order(round(as.numeric(up.gene$logFC)),decreasing = TRUE, na.last = TRUE),])
down.gene = de.gene.expr[which(de.gene.expr$logFC < 0), ]
head(down.gene[order(round(as.numeric(down.gene$logFC)),decreasing = FALSE, na.last = TRUE),])
down.gene = down.gene[order(round(as.numeric(down.gene$logFC)),decreasing = FALSE, na.last = TRUE),]
write.table(up.gene, file="up_gene_expr.csv", sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(down.gene, file="down_gene_expr.csv", sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)


library(GOstats)
all.ids <- eset.filt$ENTREZID
up.ids <- up.gene$ENTREZID
down.ids <- down.gene$ENTREZID

up.bp.params <- new(
  "GOHyperGParams", geneIds = up.ids, universeGeneIds = all.ids, 
  annotation = "hugene10sttranscriptcluster.db",
  ontology = "BP", pvalueCutoff = 0.05, conditional = FALSE,
  testDirection = "over")
up.bp.over <- hyperGTest(up.bp.params)
up.bp.over.df = summary(up.bp.over)
down.bp.params <- new(
  "GOHyperGParams", geneIds = down.ids, universeGeneIds =
    all.ids, annotation = "hugene10sttranscriptcluster.db",
  ontology = "BP", pvalueCutoff = 0.05, conditional = FALSE,
  testDirection = "over")
down.bp.over <- hyperGTest(down.bp.params)
down.bp.over.df = summary(down.bp.over)

write.table(up.bp.over.df, file="Age_up_gene.GOstats.enriched_BP_term.csv",
            row.names = FALSE, sep=",", quote = FALSE)
write.table(down.bp.over.df, file="Age_down_gene.GOstats.enriched_BP_term.csv",
            row.names = FALSE, sep=",", quote = FALSE)

g<-termGraphs(down.bp.over, id = NULL, pvalue = NULL, use.terms = TRUE)
g<-inducedTermGraph(down.bp.over,  children = TRUE, parents = TRUE)
plotGOTermGraph(g, down.bp.over = NULL, add.counts = TRUE, max.nchar = 20,
                node.colors=c(sig="lightgray", not="white"),node.shape="plaintext")

par(mfrow=c(3,3))

lapply(g, plotGOTermGraph, down.bp.over,
       node.colors=c(sig="white", not="blue"),
       node.shape="circle", add.counts=FALSE)
goplot(up.bp.over)


###################################  GO Enrichment    ##########################################


library(ALL)
library(hgu95av2.db)
library(GO.db)
library(annotate)
library(genefilter)
library(GOstats)
library(RColorBrewer)
library(xtable)
library(Rgraphviz)
library(KEGG.db)
library(org.Hs.eg.db)
library(GSEABase)
library(magrittr)
library(GOstats)
library(hgu95av2.db)
library(GO.db)
library(clusterProfiler)


library(enrichplot)



frame = toTable(org.Hs.egPATH)
keggframeData = data.frame(frame$path_id, frame$gene_id)
head(keggframeData)
keggFrame=KEGGFrame(keggframeData,organism="Homo sapiens")

data(geneList, package="DOSE")
de <- names(geneList)[abs(geneList) > 2]
de <- names(geneList)
vect<-vector()
for (i in de.gene.expr$ENTREZID) {
  if (i   %in% de   ) {
    vect<-c(vect,i)
  }else{
    print(i)
  }
  
}


#de<-de.gene.expr$ENTREZID
ego <- enrichGO(gene=as.vector(de), OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
ego <- enrichGO(gene=vect[1:5], OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)

#ego <- enrichGO(up.gene$ENTREZID[1:100], keyType = "ENTREZID", OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)

head(as.data.frame(ego))


goplot(ego, showCategory = 20)

png(file = "GoBar.png", width = 3000, height
    = 4000, res = 600)
barplot(ego, showCategory = 20)
dev.off()

enrichplot::dotplot(ego, showCategory=20)
ggplot(data = ego, aes(x = GeneRatio, y = ID, 
                       color = `p.adjust`, size = GeneRatio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ggtitle("GO enrichment analysis")
GOplot::GOBubble(ego2, labels=3)
