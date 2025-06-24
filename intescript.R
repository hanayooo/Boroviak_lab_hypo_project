rm(list = ls())
setwd("/Users/hanayoo/Desktop/cam/humansinglecell")
#install.packages("Seurat")

library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)

metatyser <- readRDS("/Users/hanayoo/Desktop/cam/single cell sequencing data/Tyser, 2021, nature/umap.rds")
write.table(metatyser, "metaTyser.txt",sep = "\t")

Guo <- readRDS("TPMGuo.RDS")
metaGuo <- read.csv("originGuo.csv", row.names = 1)
colnames(Guo) <- row.names(metaGuo)
Guototal <- CreateSeuratObject(counts = Guo, project = "Guo", min.cells = 3, min.features = 200)
Guototal <- AddMetaData(object = Guototal, metadata = metaGuo, col.name = "Origin" )
Guo <- subset(Guototal, subset = Origin == "Blastocyst")
write.table(colnames(Guo), "colguo.txt",sep = "\t")

metaGuo$cellname <- colnames(Guo)
saveRDS(Guo, "Guoseurat.RDS")

Guo <- readRDS("Guoseurat.RDS")
Liu <- read.csv("TPMLiu.csv", header = T, row.names = 1)
Petro <- read.table("TPMPetro.txt", header = T)
Tyser <- read.table("TPMTyser.txt", header = T)
Wang <- read.table("TPMWang.txt", header = T)
Xiang <- read.table("TPMXiang.txt", header = T)
Zhou <- read.table("TPMZhou.txt", header = T)

Liu <- CreateSeuratObject(counts = Liu, project = "Liu", min.cells = 3, min.features = 200)
Petro <- CreateSeuratObject(counts = Petro, project = "Petro", min.cells = 3, min.features = 200)
Tyser <- CreateSeuratObject(counts = Tyser, project = "Tyser", min.cells = 3, min.features = 200)
Wang <- CreateSeuratObject(counts = Wang, project = "Wang", min.cells = 3, min.features = 200)
Xiang <- CreateSeuratObject(counts = Xiang, project = "Xiang", min.cells = 3, min.features = 200)
Zhou <- CreateSeuratObject(counts = Zhou, project = "Zhou", min.cells = 3, min.features = 200)

saveRDS(sc.combined,"7humanplusmarmoset.RDS")
saveRDS(sc.anchors,"scanchor.RDS")
sc.list <- readRDS("sclist.RDS")
sc.list <- c(Liu, Petro, Tyser,Wang,Xiang,Zhou,Guo,marmoset)
sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "RC", scale.factor = 10000)
  
})


# select features that are repeatedly variable across datasets for integration
?SelectIntegrationFeatures
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 20000)
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, anchor.features = features)
?IntegrateData
sc.combined <- IntegrateData(anchorset = sc.anchors, k.weight = 30) 
sc.combined <- FindVariableFeatures(sc.combined, selection.method = "vst", nfeatures = 50000)
DefaultAssay(sc.combined) <- "RNA"

meta <- read.table("cells.txt", sep = "\t", header = TRUE, row.names = 1)
metaXiang <- read.table("metaXiang.txt",sep = "\t", header = TRUE, row.names = 1)

sc.combined <- AddMetaData(object = sc.combined, metadata = meta, col.name = c("author","oriano","orianos","day","seq","orianos2","subcluster") )


write.table(colnames(sc.combined@assays[["integrated"]]@data), "cells.txt",sep = "\t")
"%nin%" = Negate("%in%")
#sc.combined <- subset(sc.combined, subset = orianos %nin% c("Tb_CS3","Tb_CS5","Tb_CS6","Tb_CS7","Inner Cell Mass-Trophectoderm Transition","ISK","MIX","HE","Early Trophectoderm","Unknown","TE", "MTB","EVT","CTB","not applicable","STB","pFCC"))

#240828-resort#
sc.combined <-readRDS("CPM7humanplusmarmoset.RDS")
sc.combined <- subset(sc.combined, subset = orianos %nin% c("Gland_CS5","Gland_CS6","Gland_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","Myo_CS7","ReGland_CS5","ReGland_CS7","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Tb_CS3","Tb_CS5","Tb_CS6","Tb_CS7","Inner Cell Mass-Trophectoderm Transition","ISK","MIX","HE","Early Trophectoderm","Unknown","TE", "MTB","EVT","CTB","not applicable","STB","pFCC"))
DefaultAssay(sc.combined) <- "integrated"
saveRDS(sc.combined, file = "CPM7humanplusmamembryoonly240828.RDS")
saveRDS(sc.combined, file = "CPM7humanplusmarmoset.RDS")


sc.list <- readRDS("mouse7.5list.RDS")


sc.combined <- NormalizeData(sc.combined, normalization.method = "RC", scale.factor = 10000)
sc.combined <- ScaleData(sc.combined, verbose = FALSE)
sc.combined <- RunPCA(sc.combined, npcs = 30, verbose = FALSE)
sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindNeighbors(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindClusters(sc.combined, resolution = 0.4)
subset
?FeaturePlot

?FeaturePlot
p0 <- FeaturePlot(sc.combined, features = c("NLRP7"), reduction = "pca",min.cutoff = 0,max.cutoff = 20)
p1 <- FeaturePlot(sc.combined, features = c("EGLN3"),reduction = "pca", min.cutoff = 0,max.cutoff = 20)
p2 <- FeaturePlot(sc.combined, features = c("PA2G4"),reduction = "pca", min.cutoff = 0,max.cutoff =10)
p3 <- FeaturePlot(sc.combined, features = c("POU5F1"), reduction = "pca",min.cutoff = 0,max.cutoff = 10)
p4 <- FeaturePlot(sc.combined, features = c("CD53"),reduction = "pca", min.cutoff = 0,max.cutoff = 5)
p5 <- FeaturePlot(sc.combined, features = c("VTN"), min.cutoff = 0,max.cutoff = 2)
p6 <- FeaturePlot(sc.combined, features = c("CDKN1C"), min.cutoff = 0,max.cutoff = 4)
p7 <- FeaturePlot(sc.combined, features = c("GSTA1"), min.cutoff = 0,max.cutoff = 2)

p8 <- FeaturePlot(sc.combined, features = c("EPAS1"), min.cutoff = 0,max.cutoff = 2)
p9 <- FeaturePlot(sc.combined, features = c("FABP1"), min.cutoff = 0,max.cutoff = 2)
p10 <- FeaturePlot(sc.combined, features = c("SPINK1"), min.cutoff = 0,max.cutoff = 2)
FeaturePlot(sc.combined, features = c("FAU","EIF2S1","MAEA","PSMC1","H2AFZ"),reduction = "pca")

FeaturePlot(marmoset, features = c("AFP"),reduction = "umap",min.cutoff = 0,max.cutoff = 4)
#p12 <- FeaturePlot(sc.combined, features = c("NODAL","ESAM","FN1","FST","LHX1","TRPA1"), reduction = "pca",min.cutoff = 0,split.by = "species",ncol = 4)
#p12

VlnPlot(sc.combined, feature = "NID2", group.by = "orianos2")

saveRDS(p,"p.RDS")
saveRDS(pp,"pp.RDS")
p <- wrap_plots(p0,p1,p2,p3,p4,p5,p6,p7, ncol = 4)
p <- wrap_plots(p0,p1,p2,p3, ncol = 4)
p
p0+p1+p2
#write.table(sc.combined@reductions[["pca"]]@cell.embeddings,"PCA.txt", sep = "\t")
p0
save
?DimPlot
?WhichCells
p1 <- DimPlot(sc.combined, reduction = "umap", group.by = "day2",label = T, label.size = 2, pt.size = 0.2)
p2 <- DimPlot(sc.combined, reduction = "umap", group.by = "orianos2",label = T, label.size = 2,pt.size = 0.2, dims = c(1,2))
p3 <- DimPlot(sc.combined, reduction = "umap",  pt.size = 0.2,label = T, label.size = 4)
p1
p3

?saveRDS
save
sc.combined <- readRDS("CPM7humanplusmamembryoonly240513.RDS")
p

pdf("240513hypolinmarker.pdf", width = 12, height = 2)
plot(p)
dev.off()

