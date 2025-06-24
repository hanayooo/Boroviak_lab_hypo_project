rm(list = ls())
setwd("/Users/hanayoo/Desktop/cam/humansinglecell")

library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(ggplot2)

#embryo only#
sc.combined <- readRDS("CPM7humanplusmamembryoonly240828.RDS")
meta <- read.table("cells.txt", sep = "\t", header = TRUE, row.names = 1)
sc.combined <- AddMetaData(object = sc.combined, metadata = meta, col.name = c("author","oriano","orianos","day","seq","orianos2","subcluster","species","day2","annos") )
#rename rough idents#
#new.cluster.ids<-c("Other embryonic lineages", 
                   "Other embryonic lineages", 
                   "Hypoblast and Endoderm lineages", 
                   "Other embryonic lineages", 
                   "Other embryonic lineages", 
                   "Other embryonic lineages", 
                   "Other embryonic lineages", 
                   "Marmoset pre-implantation", 
                   "Other embryonic lineages", 
                   "Hypoblast and Endoderm lineages", 
                   "Other embryonic lineages", 
                   "Other embryonic lineages", 
                   "Marmoset pre-implantation", 
                   "Other embryonic lineages", 
                   "Other embryonic lineages")
#names(new.cluster.ids)<-levels(sc.combined)
#sc.combined<-RenameIdents(sc.combined, new.cluster.ids)
DimPlot(sc.combined)

DefaultAssay(sc.combined) <- "RNA"
Figure1A <- DimPlot(sc.combined, reduction = "umap" ,group.by = "species",pt.size = 0.2,label = T, label.size = 4)
FigureS1B <- DimPlot(sc.combined, reduction = "umap" ,group.by = "orianos2",pt.size = 0.2,label = T, label.size = 2)
Figure1B <- DimPlot(sc.combined, reduction = "umap" ,pt.size = 0.2,label = F)
Figure1C <- FeaturePlot(sc.combined, reduction = "umap",feature = c("GATA6","GATA4","PDGFRA","SOX17"), pt.size = 0.5,min.cutoff = 0,max.cutoff = 2, cols = c("grey","#0E53AB"))
A <- FeaturePlot(sc.combined, reduction = "umap",feature = c("GATA6"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 6, cols = c("grey","#0E53AB"))
A

FeaturePlot(sc.combined, reduction = "umap",feature = c("HAND1"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 1, cols = c("grey","#0E53AB"))
VlnPlot(sc.combined,features = c("SOX2","POU5F1"), group.by = "orianos2",layer = "data")
Figure1C
FeaturePlot(sc.combined, reduction = "umap",feature = c("HAND1","HAND2","FOXF1","LUM","COL3A1"), pt.size = 0.5,min.cutoff = 0,max.cutoff = 2, cols = c("grey","#0E53AB"))
a <- FeaturePlot(sc.combined, reduction = "umap",feature = c("HAND1"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 10, cols = c("grey","#0E53AB"))
b <- FeaturePlot(sc.combined, reduction = "umap",feature = c("FOXF1"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 1, cols = c("grey","#0E53AB"))
c <- FeaturePlot(sc.combined, reduction = "umap",feature = c("HAND2"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 1, cols = c("grey","#0E53AB"))
d <- FeaturePlot(sc.combined, reduction = "umap",feature = c("COL3A1"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 1, cols = c("grey","#0E53AB"))
e <- FeaturePlot(sc.combined, reduction = "umap",feature = c("LUM"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 1, cols = c("grey","#0E53AB"))
f <- FeaturePlot(sc.combined, reduction = "umap",feature = c("DCN"), pt.size = 0.2,min.cutoff = 0,max.cutoff = 10, cols = c("grey","#0E53AB"))
pp <- wrap_plots(a,c,b,f,e,d,ncol = 3)
pp
f
#select hypolin#
sc.combined <- subset(sc.combined, idents = c("Hypoblast and Endoderm lineages","Marmoset pre-implantation"))
sc.combined <- subset(sc.combined, subset = orianos2 %in% c("Hyp_CS3", "HYPO","Endoderm","SYS_CS5","SYS_CS6","SYS_CS7","VE_CS5","VE_CS6","VE_CS7","DE(NP)","DE(P)","Hypoblast","YS Endoderm"))

sc.combined <- ScaleData(sc.combined, verbose = FALSE)
sc.combined <- RunPCA(sc.combined, npcs = 30, verbose = FALSE)

sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- RunTSNE(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindNeighbors(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindClusters(sc.combined, resolution = 2)


#hypolin#
hypolin <- readRDS("/Users/hanayoo/Desktop/cam/humansinglecell/hypolin0828.RDS")
meta <- read.table("/Users/hanayoo/Desktop/cam/humansinglecell/cells.txt", sep = "\t", header = TRUE, row.names = 1)
hypolin <- AddMetaData(object = hypolin, metadata = meta, col.name = c("author","oriano","orianos","day","seq","orianos2","subcluster","species","day2","annos") )
DefaultAssay(hypolin) <- "RNA"
VlnPlot(hypolin,features =c("FN1","AFP") ,slot = "data",group.by = "day2")
FeaturePlot(hypolin, features = c("HAND1","FN1"),reduction = "pca")
DefaultAssay(hypolin) <- "RNA"

Figure1D <- DimPlot(hypolin, reduction = "pca" ,group.by = "orianos2",pt.size = 0.2,label = F)
DimPlot(hypolin, reduction = "pca" ,group.by = "orianos2",pt.size = 0.2,label = F)

FigureS1D <- DimPlot(hypolin, reduction = "pca" ,group.by = "author",pt.size = 0.2,label = F)
FigureS1C <- DimPlot(hypolin, reduction = "pca" ,group.by = "day2",pt.size = 0.2,label = F)
  Figure1E1 <- FeaturePlot(hypolin, reduction = "pca",feature = c("EOMES"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 2,cols = c("grey","#0E53AB"),label = FALSE,)
  Figure1E2 <- FeaturePlot(hypolin, reduction = "pca",feature = c("LEFTY2"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 5,cols = c("grey","#0E53AB"))
  Figure1E3 <- FeaturePlot(hypolin, reduction = "pca",feature = c("NODAL"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 5,cols = c("grey","#0E53AB"))
  Figure1E4 <- FeaturePlot(hypolin, reduction = "pca",feature = c("TTR"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 20,cols = c("grey","#0E53AB"))
  Figure1E5 <- FeaturePlot(hypolin, reduction = "pca",feature = c("APOB"), pt.size = 0.5,min.cutoff = 5, max.cutoff = 40,cols = c("grey","#0E53AB"))
  Figure1E6 <- FeaturePlot(hypolin, reduction = "pca",feature = c("RBP4"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 2,cols = c("grey","#0E53AB"))
Figure1E <- wrap_plots(Figure1E1,Figure1E2,Figure1E3,Figure1E4,Figure1E5,Figure1E6,ncol = 3)

DefaultAssay(hypolin) <- "RNA"

###
hypolin <- FindClusters(hypolin, resolution = 2)
new.cluster.ids<-c("Endo", "Pre_HYPO", "Endo", "Pre_HYPO", "SYSE","Pre_HYPO","VE","Endo","VE","Endo","SYSE"
)
names(new.cluster.ids)<-levels(hypolin)
hypolin<-RenameIdents(hypolin, new.cluster.ids)
Figure1F <- DimPlot(hypolin, reduction = "pca",pt.size = 0.2,label = F)
Figure1F

sc.combinedhuman <- subset(hypolin, subset = species == "human")
DefaultAssay(sc.combinedhuman) <- "RNA"

markers <- FindAllMarkers(sc.combinedhuman,logfc.threshold = 1,min.pct = 0.5,only.pos = T)
markers_filtered <- subset(markers, markers$p_val_adj < 0.05)
VEmarkers <- subset(markers_filtered, markers_filtered$cluster == "VE")
VEmarkers$ratio <- VEmarkers$pct.1/VEmarkers$pct.2
VEmarkerst100 <- VEmarkers %>% top_n(100,avg_log2FC)

VE <- subset(hypolin, idents = "VE" )
humanVE <- subset(VE, subset = species == "human" )

humanVEexpression <- humanVE@assays[["RNA"]]@data
humanVEexpressiont <-t(humanVEexpression)
humanVEexpressiont <- as.data.frame(humanVEexpressiont)

humanVEgenes <- humanVEexpressiont[VEmarkers$gene]
cor (humanVEgenes, method="pearson")
humantdc <- cor (humanVEgenes, method="pearson")
corrplot(humantdc)


# 提取多个 AVE marker 的相关性
ave_genes <- c("CER1", "LEFTY1", "LEFTY2")
ave_corr_matrix <- humantdc[ave_genes, ]

# 计算每个基因与这3个基因的平均相关性
mean_corr <- colMeans(ave_corr_matrix, na.rm = TRUE)

# 排序
sorted_by_ave <- sort(mean_corr, decreasing = TRUE)
gene_order <- names(sorted_by_ave)
ordered_cor_matrix <- humantdc[gene_order, gene_order]
corrplot(ordered_cor_matrix) 

# 查看结果
head(sorted_by_ave, 20)

# 假设你的相关性矩阵是 cor_matrix
dist_mat <- as.dist(1 - ordered_cor_matrix)        # 距离矩阵，相关性越高距离越近
hc <- hclust(dist_mat, method = "ward.D2") # 层次聚类

# 切成 4 个 cluster（你可以改成 5、6 等）
clusters <- cutree(hc, k = 6)

# 把 cluster 信息做成注释
annotation <- data.frame(Cluster = as.factor(clusters))
rownames(annotation) <- names(clusters)

# 用 pheatmap 画图，加上 annotation
pheatmap(ordered_cor_matrix,
         cluster_rows = hc,
         cluster_cols = hc,
         annotation_row = annotation,
         annotation_col = annotation,
         show_rownames = T,
         fontsize = 3,
         show_colnames = T)
ordered_cor_matrix

library(pheatmap)
pheatmap(ordered_cor_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = T,
         show_colnames = T,
         fontsize = 3,
         main = "Correlation heatmap ")
?pheatmap


###heatmap
PreHYPO10 <- read.table("prehypo20.txt")
VEmarker <- c("CER1","LEFTY2","HHEX","LHX1","FZD5","NODAL","FN1","LGR5","HAS2","TDGF1","CXCR4","COL4A5","FZD7","OTX2","LRIG3","EOMES","SHISA2","LRRN1","GPC4")
VEmarker10 <- c("CER1","LEFTY2","HHEX","LHX1","FZD5","NODAL","FN1","EOMES","OTX2","CXCR4")
SYSmarker <-c("TTR","TF","APOB","RBP4","SERPINA1","VTN","CDKN1C","GSTA1","EPAS1","FABP1","SPINK1")
SYSmarker10 <-c("TTR","TF","APOB","RBP4","SERPINA1","VTN","CDKN1C","GSTA1","EPAS1","FABP1")


?FindVariableFeatures
hypolin <- FindVariableFeatures(hypolin, nfeatures = 50000)
hypolin <- ScaleData(hypolin)
hypolinheatmap <- subset(hypolin, idents = c("Pre_HYPO","SYSE","VE"))
DefaultAssay(hypolinheatmap) <- "RNA"
P <- DoHeatmap(hypolinheatmap,features = c(PreHYPO10$V1,SYSmarker10,VEmarker10),slot = "scale.data",size = 4)+ scale_fill_gradientn(colors = c("#0E53AB", "#eeeeee", "#ef4242"))
P<-corrplot(tdc)
P
??DoHeatmap

pdf("250412VEmarkers.pdf", width = 12, height = 12)
pheatmap(ordered_cor_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = T,
         show_colnames = T,
         fontsize = 3,
         main = "Correlation heatmap ")
dev.off()
FigureS
#DoHeatmap(CER1,features = marker.list$gene,slot = "scale.data",size = 4)+ scale_fill_gradientn(colors = c("#0E53AB", "#eeeeee", "#ef4242"))
total
sc.combined <-readRDS("CPM7humanplusmarmoset.RDS")
total <- DimPlot(sc.combined, cells.highlight = list("TE" = TE,"Maternal"= Maternal,"other"= other),cols.highlight = c("#67A1C5","#0E53AB","#F6B46D"))
?DimPlot
?WhichCells
TE <- rownames(meta[meta$oriano %in% c("Tb_CS3","Tb_CS5","Tb_CS6","Tb_CS7","Inner Cell Mass-Trophectoderm Transition","Early Trophectoderm","TE", "MTB","EVT","CTB","STB","pFCC"),])
other <- rownames(meta[meta$oriano %in% c("not applicable","Unknown","ISK","MIX","HE"),])
Maternal <- rownames(meta[meta$oriano %in% c("Gland_CS5","Gland_CS6","Gland_CS7","Stroma_CS5","Stroma_CS6","Stroma_CS7","Myo_CS7","ReGland_CS5","ReGland_CS7","ReStroma_CS5","ReStroma_CS6","ReStroma_CS7"),])


hypolin <- RunPCA(hypolin,npcs = 30)

pca.stdev <- hypolin@reductions$pca@stdev  # 各PC的标准差
pca.var <- pca.stdev^2                 # 方差 = stdev^2
var.percent <- pca.var / sum(pca.var) * 100  # 计算每个PC的方差百分比

# 查看前10个PC的方差解释率
head(var.percent, 10)
ElbowPlot(hypolin)

FeaturePlot(hypolin, reduction = "pca",feature = c("APOA4"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 1,cols = c("grey","#0E53AB"),label = FALSE,)
FeaturePlot(hypolin, reduction = "pca",feature = c("VTN"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 1,cols = c("grey","#0E53AB"),label = FALSE,)
FeaturePlot(hypolin, reduction = "pca",feature = c("TTR"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 1,cols = c("grey","#0E53AB"),label = FALSE,)
FeaturePlot(hypolin, reduction = "pca",feature = c("APOB"), pt.size = 0.5,min.cutoff = 0, max.cutoff = 20,cols = c("grey","#0E53AB"),label = FALSE,)


pdf("250616exmesfeatureplot.pdf", width = 12, height = 7)
plot(pp)
dev.off()
