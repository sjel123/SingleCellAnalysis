library(Seurat)
library(Matrix)
library(dplyr)
library(RColorBrewer)
library(colorRamps)
# Load in dataset 
#setwd("~/projects/Colonna_SC/")
#Load Seurat
load ("../Data/Lupus.data.Robj")
Metadata <- read.table("../Data//metadata.csv", header=TRUE, sep=",")
Data <- read.table("../Data/METRO_800_Post-QC.csv", header =TRUE, sep=",")


#Filter Data to limit those genes with a lot of missing values
#data5Fil <- data5[which(apply(data5[,2:length(data5)], 1, function(x) sum(x >0))>1400),]
#data5Fil <- data5[apply(data5[,2:length(data5)], 1, function(x) sum(x)>1300),]

# transform data to log scale
#data1.log=log(data1[,2:382]+1)

#Create spearse matrix to load into Seurat frame
m = Matrix(as.matrix(Data[,2:length(Data)]))
  m <- as(m, "dgTMatrix")
  row.names(m) <- Data[,1]



Metadata.f <- Metadata[Metadata$SAMPLE_ID %in% colnames(m),]
table(Metadata.f$SAMPLE_ID==colnames(m))


# m = Matrix(as.matrix(data2[,-1]))
#   m <- as(m, "dgTMatrix")
#   row.names(m) <- data2[,1]
# Initialize the Seurat object with the raw (non-normalized data)
Lupus.data <- new("seurat", raw.data=m)

# Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
# Perform log-normalization, first scaling each cell to a total of 1e4 molecules 
# (as in Macosko et al. Cell 2015)
Lupus.data <- Setup(Lupus.data, min.cells = 3, min.genes = 200, do.logNormalize =T, total.expr = 1e4, project = "AMP_Lupus")


#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
Samples <- data.frame(Lupus.data@data@Dimnames[2])
Metadata.f1 <- merge(Samples, Metadata, by.x=1, by.y="SAMPLE_ID")
  Metadata.f1<- Metadata.f1[match(rownames(Lupus.data@data.info),Metadata.f1[,1]),]
  table(Metadata.f1[1]==row.names(Lupus.data@data.info))

  X <- Metadata.f1$X
  attr(X, "names")<- as.character(Metadata.f1[,1] )
Subject <- Metadata.f1$SUBJECT_ID
  attr(Subject, "names")<- as.character(Metadata.f1[,1] )
Tissue <- Metadata.f1$TISSUE
  attr(Tissue, "names")<-as.character(Metadata.f1[,1])


Lupus.data <- AddMetaData(Lupus.data, X, "X")
Lupus.data <- AddMetaData(Lupus.data, Subject, "Subject")
Lupus.data <- AddMetaData(Lupus.data, Tissue, "Tissue")

#Clean up
rm(Data)
rm(list=ls(pattern="Metadat*"))
rm(m)
rm(Samples)

#nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, 
#nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito 
#using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
mito.genes <- grep("^MT-", rownames(Lupus.data@data), value = T)
percent.mito <- colSums(expm1(Lupus.data@data[mito.genes, ]))/colSums(expm1(Lupus.data@data))

  
#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
  Lupus.data <- AddMetaData(Lupus.data, percent.mito, "percent.mito")

  #Lupus.data <- AddMetaData(Lupus.data, plate, "Plate")
  VlnPlot(Lupus.data, c("nGene", "nUMI", "percent.mito"), nCol = 3)

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(Lupus.data, "nUMI", "percent.mito")
GenePlot(Lupus.data, "nUMI", "nGene")

#We filter out cells that have unique gene counts over 2,500 and under 500, and > 5% mitochondrial percentage
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)

#Lupus.data <- SubsetData(Lupus.data, subset.name = "nGene", accept.high = 5000)
#Lupus.data <- SubsetData(Lupus.data, subset.name = "percent.mito", accept.high = 0.8)
Lupus.data <- SubsetData(Lupus.data , subset.name = "nGene", accept.low = 600)

#choose gene outliers on mean-variability plot
#Detection of variable genes across the single cells
Lupus.data <- MeanVarPlot(Lupus.data, x.low.cutoff = 0.1, y.cutoff = 0.5)
#Lupus.data <- MeanVarPlot(Lupus.data, x.low.cutoff = 0.1, y.cutoff = 1)
  length(Lupus.data@var.genes)
  #Lupus.data <- MeanVarPlot(Lupus.data, fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
    #length(Lupus.data@var.genes)

#Perform negative-binomial regression on the variable genes, this sets their value in 
#Lupus.data@scale.data, which is used for PCA/clustering
#We only do this on the variable genes to save time, but you can do this genome-wide
#We treat mitochondrial percentage, batch, and nUMI as confounding variables, 

#you can save the object at any time to save results, and can restore it back in using load()
  save(Lupus.data, file = "../Data/Lupus.data.Robj")

#Regress out unwanted sources of variation
Lupus.data <- RegressOut(Lupus.data, latent.vars = c ("nUMI"))#, "percent.mito"))
#Lupus.data <- RegressOut(Lupus.data, latent.vars = c("percent.mito", "orig.ident", "nUMI"), genes.regress = Lupus.data@var.genes, model.use = "negbinom")

#Run PCA with the IRLBA package (iteratively computes the top dimensions, dramatic increase in speed since we are throwing away most PCs anyway)
#Lupus.data <- PCAFast(Lupus.data, pc.genes = Lupus.data@var.genes, pcs.compute = 40, pcs.print = 30)
Lupus.data <- PCA(Lupus.data, pc.genes = Lupus.data@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)

PCElbowPlot(Lupus.data, num.pc = 40)
  PrintPCA(Lupus.data,pcs.print = 1:36)
  PCHeatmap(Lupus.data,pc.use = 1:12)
  PCHeatmap(Lupus.data,pc.use = 13:24,100)
  PCHeatmap(Lupus.data,pc.use = 25:36,100)

#Lupus.data@data.info <- merge(Lupus.data@data.info[,c(1:2,7:8)], Metadata, by.x=0, by.y="SAMPLE_ID", all.x=TRUE)


#select 3 PCs for downstream analysis

for (i in seq(2,40, by=3)){print(i)
  Lupus.data <- RunTSNE(Lupus.data, dims.use = 1:i, do.fast = T)
  assign(paste0("TSNPlot",i), TSNEPlot(Lupus.data,do.label = T, group.by="Tissue",do.return=T))
}
for (i in seq(2,40, by=3)){
 print(paste0("TSNPlot",i,"$labels$title=TSNE",i))
}

for (i in seq(2,40, by=3)){
  p = get(paste0("TSNPlot",i))
  p+ggtitle(paste0("TSNE",i))
  p
}

TSNPlot2$labels$title="TSNE2"
TSNPlot5$labels$title="TSNE5"
TSNPlot8$labels$title="TSNE8"
TSNPlot11$labels$title="TSNE11"
TSNPlot14$labels$title="TSNE14"
TSNPlot17$labels$title="TSNE17"
TSNPlot20$labels$title="TSNE20"
TSNPlot23$labels$title="TSNE23"
TSNPlot26$labels$title="TSNE26"
TSNPlot29$labels$title="TSNE29" 
TSNPlot32$labels$title="TSNE32"
TSNPlot35$labels$title="TSNE35"
TSNPlot38$labels$title="TSNE38"
  
for (i in seq(2,40, by=3)){
  print(i)
  print(get(paste0("TSNPlot",i)))
}

write.table(TSNPlot5$data, "../Data/TNSE.txt", sep="\t")
write.table(as.data.frame(as.matrix(Lupus.data@data)), "../Data/ExprsData.txt", sep="\t")



i=5
Lupus.data <- RunTSNE(Lupus.data, dims.use = 1:i, do.fast = T)

#save.SNN means that you can easily re-run with different resolution values. Here we run with a few different res v
#Lupus.data <- FindClusters(Lupus.data ,pc.use = 1:5, resolution = 0.4, save.SNN = T, do.sparse = T, temp.file.location = "./")
#Lupus.data <- FindClusters(Lupus.data ,pc.use = 1:38, resolution = 1.2, save.SNN = T, do.sparse = T, temp.file.location = "./")
Lupus.data <- FindClusters(Lupus.data ,pc.use = 1:5, resolution = 0.2, save.SNN = T, do.sparse = T, temp.file.location = "./")

#TSNEPlot(Lupus.data,do.label = T, group.by="res.0.4")
TSNEPlot(Lupus.data,do.label = T, group.by="res.0.2")
TSNEPlot(Lupus.data,do.label = T, group.by="Tissue")
TSNEPlot(Lupus.data,do.label = T, group.by="Subject")
TSNEPlot(Lupus.data,do.label = T, group.by="Subject")
TSNEPlot(Lupus.data,do.label = T, group.by="X")

# Lupus.data <- FindClusters(Lupus.data ,pc.use = 1:2, resolution = 0.2, save.SNN = T, do.sparse = T, temp.file.location = "./")
# TSNEPlot(Lupus.data,do.label = T, group.by="res.0.2")
# Lupus.data <- FindClusters(Lupus.data ,pc.use = 1:5, resolution = 0.8, save.SNN = T, do.sparse = T, temp.file.location = "./")

#Lupus.data_20=FindClusters(Lupus.data,pc.use = 1:12,resolution = seq(0.6,4,0.1),save.SNN = T,do.sparse = T,k.param = 20, temp.file.location = "./")
#Lupus.data <- FindClusters(Lupus.data ,pc.use = 1:3, resolution = seq(0.3,1.5,0.1), save.SNN = T, do.sparse = T, temp.file.location = "./")
save(Lupus.data, file = "Data/Lupus.data1.Robj")

#Switch labels if multiple clusters identified
#Search at Lupus.data@data.info$
Lupus.data <- SetAllIdent(Lupus.data, id = "res.0.4")
#Lupus.data <- SetAllIdent(Lupus.data, id = "Plate")
TSNEPlot(Lupus.data,do.label = T, group.by="res.0.6")
Lupus.data <- SetAllIdent(Lupus.data,id = "TISSUE")
Lupus.data <- SetAllIdent(Lupus.data,id = "res.0.2")
TSNEPlot(Lupus.data,do.label = T)

#Lupus.data <- PCA(Lupus.data, pc.genes = Lupus.data@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
#ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation with the calculated components. Though we don’t use this further here, it can be used to identify markers that are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. The results of the projected PCA can be explored by setting use.full=T in the functions below.

Lupus.data <- ProjectPCA(Lupus.data)

#Seurat provides several useful ways of visualizing both cells and genes that define the PCA, including PrintPCA(), VizPCA(), PCAPlot(), and PCHeatmap()
par(mfrow = c(1, 1))
# Examine  and visualize PCA results a few different ways
PrintPCA(Lupus.data, pcs.print = 1:12, genes.print = 5, use.full = TRUE)
VizPCA(Lupus.data, 1:4)
PCAPlot(Lupus.data, 1, 2)
PCAPlot(Lupus.data, 1, 3)
PCHeatmap(Lupus.data, pc.use = 1, cells.use = 100, do.balanced = TRUE)
PCHeatmap(Lupus.data, pc.use = 1:12, cells.use = 100, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(Lupus.data, pc.use = 13:24, cells.use = 100, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)


# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
Lupus.data <- JackStraw(Lupus.data, num.replicate = 100, do.print = FALSE)
JackStrawPlot(Lupus.data, PCs = 1:12)
PCElbowPlot(Lupus.data, num.pc = 20)


#Run Non-linear dimensional reduction (tSNE)

#Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets. While we no longer advise clustering directly on tSNE components, cells within the graph-based clusters determined above should co-localize on the tSNE plot. This is because the tSNE aims to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. As input to the tSNE, we suggest using the same PCs as input to the clustering analysis, although computing the tSNE based on scaled gene expression is also supported using the genes.use argument.

# Lupus.data <- RunTSNE(Lupus.data, dims.use =1:2, do.fast = T)
# 
# #note that you can set do.label=T to help label individual clusters
# TSNEPlot(Lupus.data, do.label=T)
# 
# plot(Lupus.data@tsne.rot)



# find all markers of cluster 1
cluster1.markers <- FindMarkers(Lupus.data, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(Lupus.data, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(Lupus.data, ident.1 = 0, min.pct = 0.25)
print(head(cluster1.markers, 5))
print(head(cluster2.markers, 5))
print(head(cluster3.markers, 5))

Lupus.markers <- FindAllMarkers(Lupus.data, only.pos=FALSE, min.pct =0.25, thresh.use =0.25)
top2 <- Lupus.markers %>% group_by(cluster) %>% top_n(2, avg_diff)
top10 <- Lupus.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
ILC1 <- c("GZMA", "CCL5","IFNG", "KLRD1", "NKG7", "CXCR3")
ILC2 <- c("GATA3", "RORA", "S100A6")
ILC3 <- c( "AHR", "LIF", "IL23R", 'CD300LF')
NK <- c("EOMES", "CD300A", "KLRF1", "GZMB")
ILC1nonNK <- c( "CD27","CD28", "IL6R", "TNFRSF10A")
genes <- c("CCL5", "GZMA", "TCF7", "IL10RA", "KLRD1", "CD2", "ITGAM", "LAIR1"  )
genes2 <- c("NKG7", "GNLY", "GZMA", "LST1")
genes3 <- c("GZMK", "CMC1", "LST1", "CSF2", "LYST", "CCL4")
#ILC1, 2, 3 markers
genes4 <- c("CCL5", "GZMA", "IL10RA", "FCRL3", "IL7R", "LST1")
genes5 <- c("NKG7", "GZMA", "LST1", "CSF2")
genes6 <- c(   "RORA", "GATA3")
VlnPlot(Lupus.data, top2$gene)
VlnPlot(Lupus.data, c(ILC1, ILC2, ILC3))
VlnPlot(Lupus.data, ILC1)
VlnPlot(Lupus.data, NK)
VlnPlot(Lupus.data,ILC1nonNK)
VlnPlot(Lupus.data, c(ILC1[1:2], ILC2,ILC3))
VlnPlot(Lupus.data, TopMarkers )
VlnPlot(Lupus.data, "KIT" )

TopMarkers <- c("NKG7", "GZMA", "KRT86", "LST1", "FOS", "JUN", "GATA3", "LGALS1", "HPGDS", "SOX4",
                "CD9", "CDC37L1")

FeaturePlot(Lupus.data, top2$gene, cols.use=c("yellow", "red"))
FeaturePlot(Lupus.data, TopMarkers, cols.use=c("yellow", "red"))
FeaturePlot(Lupus.data, top2$gene[1:5], cols.use=colorRampPalette(brewer.pal(11,"Spectral"))(100))
FeaturePlot(Lupus.data, top2$gene, cols.use=blue2green2red(100))
FeaturePlot(Lupus.data, ILC1, cols.use=colorRampPalette(c("grey","yellow", "red"))(50))
FeaturePlot(Lupus.data,ILC2, cols.use=c("grey", "blue"))
FeaturePlot(Lupus.data,ILC3, cols.use=c("grey", "blue"))
FeaturePlot(Lupus.data, "NKG7", cols.use=colorRampPalette(c("grey","yellow", "red"))(50))
FeaturePlot(Lupus.data, "KLRD1", cols.use=colorRampPalette(c("grey","yellow", "red"))(50))

FeaturePlot(Lupus.data,top2$gene, cols.use=(blue2green2red(100)))#cols.use=c("grey", "blue"))
FeaturePlot(Lupus.data,"NKG7", cols.use=c("grey", "blue"))
FeaturePlot(Lupus.data,"CD300LF", cols.use=c("grey", "blue"))
FeaturePlot(Lupus.data, "CD300LF", cols.use=c("yellow", "red"))
Lupus.data@raw.data@Dimnames[[1]][grep("ROR", Lupus.data@raw.data@Dimnames[[1]])]

pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10
pbmc.markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(Lupus.data, genes.use = top2$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(Lupus.data, genes.use = top20$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(Lupus.data, genes.use = top10$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)

DoHeatmap(Lupus.data, genes.use = c("CCL5", "GZMA"), order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)


#Summarize expression by clusters

Lupus.data@data.info$Labes <- Lupus.data@data.info$res.0.6
  Lupus.data@data.info$Labes <- gsub(2, "ILC2a",Lupus.data@data.info$Labes )
  Lupus.data@data.info$Labes <- gsub(0, "ILC2b",Lupus.data@data.info$Labes )
  Lupus.data@data.info$Labes <- gsub(3, "ILC2c",Lupus.data@data.info$Labes )
  Lupus.data@data.info$Labes <- gsub(1, "ILC3",Lupus.data@data.info$Labes )
  Lupus.data@data.info$Labes <- gsub(4, "ILC1",Lupus.data@data.info$Labes )
  Lupus.data@data.info$Labes <- gsub(5, "ILC?",Lupus.data@data.info$Labes )
Lupus.data <- SetAllIdent(Lupus.data,id = "Labes")
# What Proportion are in each cluster
prop.table(table(Lupus.data@ident))

cluster.averages <- AverageExpression(Lupus.data)
head(cluster.averages)


# Return this information as a Seurat object (enables downstream plotting
# and analysis)
cluster.averages <- AverageExpression(object = Lupus.data, return.seurat = TRUE)#, show.progress = TRUE)
cluster.averages@data$Diff <- apply(cluster.averages@data, 1, function(x) (max(x)/min(x)))
cluster.averages.Table <- data.frame(as.data.frame(as.matrix(cluster.averages@data)),
                            FC = apply(cluster.averages@data, 1, function(x) (max(log(x))-min(log(x+0.1)))))

cluster.averages.Table[which(cluster.averages.Table$FC==max(cluster.averages.Table$FC)),]
  topAveMark <- head(cluster.averages.Table[order(cluster.averages.Table$FC,decreasing = TRUE),],40)
      round(topAveMark,1)
cluster.averages.Table.F <- cluster.averages.Table[cluster.averages.Table$FC >1,]
head(cluster.averages.Table.F[order(cluster.averages.Table.F$ILC2c,decreasing = TRUE),],20)
round(topAveMark,1)
DotPlot(Lupus.data,c(row.names(topAveMark)[1:10]),thresh.col =1, cols.use = myPalette(low="lightgrey",high = "blue"),
        cex.use = 2, use.imputed=TRUE)
boxplot(topAveMark[1:3,1:6])
, aes())+geom_point()

# How can I plot the average expression of NK cells vs. CD8 T cells?  Pass
# do.hover = T for an interactive plot to identify gene outliers
CellPlot(object = cluster.averages, cell1 = "ILC1", cell2 = "ILC3")

# You can also plot heatmaps of these 'in silico' bulk datasets to visualize
# agreement between replicates
DoHeatmap(object = cluster.averages, genes.use = PCTopGenes(object = Lupus.data, pc.use = 2, 
                                                            do.balanced = TRUE), group.label.rot = TRUE, group.cex = 0)

#Visualize canonical markers
FeaturePlot(Lupus.data, c("MS4A1", "GNLY","CD3E","CD8A","LYZ"),cols.use = c("lightgrey","blue"),nCol = 3)
FeaturePlot(Lupus.data, c("CD27","CD3E","CD8A","LYZ"),cols.use = c("lightgrey","blue"),nCol = 3)
FeaturePlot(Lupus.data, c("GZMA","CCL5","LST1","S100A4"),cols.use = c("lightgrey","blue"),nCol = 3)
FeaturePlot(Lupus.data, c("IFNG","IL22","TIGIT"),cols.use = c("lightgrey","blue"),nCol = 3)
FeaturePlot(Lupus.data, c("IFNG","IL22","TIGIT"),cols.use = c("lightgrey","blue"),nCol = 3)

#Dot plot visualization
DotPlot(Lupus.data,c(top2$gene),cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)

##Save files for Shiny App
tsne_pca.t4.1 <- Lupus.data@ tsne.rot
PCAData.1 <-  Lupus.data@pca.rot
outliers.1 <- colnames(Lupus.data@raw.data)[!colnames(Lupus.data@raw.data)%in%colnames(Lupus.data@data)]
res.1 <- data2
save( res.1, outliers.1, PCAData.1, tsne_pca.t4.1, file ="Data/PTData.1.RDa")

maxRow.1 <- apply(res[,-1], 1, max)
maxCount.1 <- apply(res[,-1],1, function(x)sum(x>0))
Var.1 <- apply(res[,-1],1, var)

class3.1 <- Lupus.data@data.info$res.0.4
save( class3.1, maxRow.1, maxCount.1, Var.1, file ="Data/StatData.1.Rda")


##############
## Run analysis on genes expressed in at least 50% of cells
##
library(Rtsne)
data1Fil <- data1[which(apply(data1[,2:382], 1, function(x) sum(x >0))>191),]
data1Filpca <- prcomp(data1Fil[,2:382], center=TRUE, scale.=TRUE)
plot(data1Filpca$x[,1], data1Filpca$x[,2])

rtsne_out <- Rtsne(data1Fil[,2:382])
plot(rtsne_out$Y, t='n', main="BarnesHutSNE")
text(rtsne_out$Y, labels=rownames(mydata))

####Plot genes of interest Druggable Pathways from Marco 9/11/17
Cytokine_receptors <- c("IL23R", "IL4R", "IFNGR1", "IFNGR2", "IL21R", "KIT")
Cell_Surface <- c("FCRL3", "TNFRSF25", "TNFRSF21", "TNFRSF18", "LTA",  "CD226", "TIGIT", "CD96")
Enzymes_signaling <- c("JAK3", "ENTPD1","NT5E", "P2RX7", "INPP4B")
Others <- c("CD300LF", "CD160", "IKZF3")


Lupus.data@data.info$ILC0.4 <- Lupus.data@data.info$res.0.4
  Lupus.data@data.info$ILC0.4  <- gsub("2", "ILC3",Lupus.data@data.info$ILC0.4  )
  Lupus.data@data.info$ILC0.4  <- gsub("0", "ILC2",Lupus.data@data.info$ILC0.4  )
  Lupus.data@data.info$ILC0.4  <- gsub("1", "ILC1",Lupus.data@data.info$ILC0.4  )
Lupus.data <- SetAllIdent(Lupus.data,id = "ILC0.4")

#Visualize canonical markers
FeaturePlot(Lupus.data, Cytokine_receptors,cols.use = c("lightgrey","blue"),nCol = 2)
FeaturePlot(Lupus.data, Cell_Surface,cols.use = c("lightgrey","blue"),nCol = 3)
FeaturePlot(Lupus.data, Enzymes_signaling ,cols.use = c("lightgrey","blue"),nCol = 2)
FeaturePlot(Lupus.data, Others,cols.use = c("lightgrey","blue"),nCol = 2)

#Dot plot visualization
DotPlot(Lupus.data,Cytokine_receptors,cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)
DotPlot(Lupus.data,Cell_Surface,cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)
DotPlot(Lupus.data,Enzymes_signaling,cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)
DotPlot(Lupus.data,Others,cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)


VlnPlot(Lupus.data, Cytokine_receptors )
VlnPlot(Lupus.data, Cell_Surface)
VlnPlot(Lupus.data, Enzymes_signaling)
VlnPlot(Lupus.data, Others)
