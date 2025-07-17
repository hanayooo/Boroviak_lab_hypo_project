rm(list = ls())
setwd("/Users/hanayoo/Desktop/cam/humansinglecell")

library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)

Guo <- read.csv("countGuo.csv", row.names = 1)
Guo2 <- read.csv("countGuo.csv", row.names = 1)
metaGuo <- read.csv("originGuo.csv", row.names = 1)
colnames(Guo) <- row.names(metaGuo)
Guototal <- CreateSeuratObject(counts = Guo, project = "Guo", min.cells = 3, min.features = 200)
Guototal <- AddMetaData(object = Guototal, metadata = metaGuo, col.name = "Origin" )
Guo <- subset(Guototal, subset = Origin == "Blastocyst")
saveRDS(Guo, "Guocountseurat.RDS")

Liu <- read.table("CPMLiu.txt", header = T, row.names = 1)
Petro <- read.table("countPetro.txt", header = T)
Tyser <- read.table("CPMTyser.txt", header = T)
Wang <- read.csv("countWang.csv", header = T, row.names = 1)
Xiang <- read.table("countXiangrmdup.txt", header = T)
Zhou <- read.table("CPMZhou.txt", header = T)

Liu <- CreateSeuratObject(counts = Liu, project = "Liu", min.cells = 3, min.features = 200)
Petro <- CreateSeuratObject(counts = Petro, project = "Petro", min.cells = 3, min.features = 200)
Tyser <- CreateSeuratObject(counts = Tyser, project = "Tyser", min.cells = 3, min.features = 200)
Wang <- CreateSeuratObject(counts = Wang, project = "Wang", min.cells = 3, min.features = 200)
Xiang <- CreateSeuratObject(counts = Xiang, project = "Xiang", min.cells = 3, min.features = 200)
Zhou <- CreateSeuratObject(counts = Zhou, project = "Zhou", min.cells = 3, min.features = 200)
Guo <- readRDS("Guocountseurat.RDS")
marmoset <- readRDS("marmosetseurat.RDS")

write.table(rownames(marmoset@assays[["RNA"]]@counts),"marmosetgene.txt")


sc.list <- c(Liu, Petro, Tyser,Wang,Xiang,Zhou,Guo,marmoset)
sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "RC", scale.factor = 10000)
  
})

features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 20000)
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, anchor.features = features)
sc.combined <- IntegrateData(anchorset = sc.anchors)
sc.combined <- FindVariableFeatures(sc.combined, selection.method = "vst", nfeatures = 20000)

DefaultAssay(sc.combined) <- "integrated"
meta <- read.table("cells.txt", sep = "\t", header = TRUE, row.names = 1)
metaXiang <- read.table("metaXiang.txt",sep = "\t", header = TRUE, row.names = 1)

sc.combined <- AddMetaData(object = sc.combined, metadata = meta, col.name = c("author","oriano","orianos","day","seq","orianos2","subcluster") )


write.table(colnames(sc.combined@assays[["integrated"]]@data), "cells.txt",sep = "\t")
"%nin%" = Negate("%in%")

sc.combined <-readRDS("CPM7humanplusmarmoset.RDS")
sc.combined <- subset(sc.combined, subset = orianos %nin% c("Gland_CS5","Gland_CS6","Gland_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","Myo_CS7","ReGland_CS5","ReGland_CS7","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7","Tb_CS3","Tb_CS5","Tb_CS6","Tb_CS7","Inner Cell Mass-Trophectoderm Transition","ISK","MIX","HE","Early Trophectoderm","Unknown","TE", "MTB","EVT","CTB","not applicable","STB","pFCC"))
DefaultAssay(sc.combined) <- "integrated"
saveRDS(sc.combined, file = "CPM7humanplusmamembryoonly240828.RDS")
saveRDS(sc.combined, file = "CPM7humanplusmarmoset.RDS")


sc.combined <- NormalizeData(sc.combined, normalization.method = "RC", scale.factor = 10000)
sc.combined <- ScaleData(sc.combined, verbose = FALSE)
sc.combined <- RunPCA(sc.combined, npcs = 30, verbose = FALSE)
sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindNeighbors(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindClusters(sc.combined, resolution = 0.4)