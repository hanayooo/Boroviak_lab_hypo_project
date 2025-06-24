rm(list = ls())
setwd("/Users/hanayoo/Desktop/cam/seq2411-mt")

library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(ggplot2)

ss22411count <- read.csv("/Users/hanayoo/Desktop/cam/seq2411-mt/250505-ss22411countrmdup-mt.csv",row.names = 1,check.names = F)
######

ss22411 <- CreateSeuratObject(counts = ss22411count, project = "ss22411", min.cells = 3, min.features = 1000)

ss22411 <- FindVariableFeatures(ss22411, selection.method = "vst", nfeatures = 60000)

ss22411 <- NormalizeData(ss22411, normalization.method = "RC", scale.factor = 10000)
ss22411 <- ScaleData(ss22411, verbose = FALSE)
ss22411 <- RunPCA(ss22411, npcs = 30, verbose = FALSE)
ss22411 <- RunUMAP(ss22411, reduction = "pca", dims = 1:15)
ss22411 <- FindNeighbors(ss22411, reduction = "pca", dims = 1:15)
ss22411 <- FindClusters(ss22411, resolution = 0.4)

sample <- colnames(ss22411count) 
groups <- sub("_\\d+$", "", sample)
df <- data.frame(cells = sample, group = groups, stringsAsFactors = FALSE)
####

ss22411 <- AddMetaData(object = ss22411, metadata = df, col.name = c("group"))

CER1 <- subset(ss22411, subset = group %in% c("cRCER1_ACL_CER1+","CER1_primed_DE_72h","cRCER1_DE_72h","CER1_primed_ACL_120h"))
CER12 <- subset(ss22411, subset = group %in% c("cRCER1_ACL_CER1+","CER1_primed_DE_72h","cRCER1_DE_72h","CER1_primed_ACL_120h","CER1_primed_ACL_72h"))
CER1 <- ScaleData(CER1)
Idents(CER1) <- CER1@meta.data$group
CER1@meta.data$group <- factor(
  CER1@meta.data$group,
  levels = c("cRCER1_ACL_CER1+", "cRCER1_DE_72h", "CER1_primed_ACL_120h", "CER1_primed_DE_72h")  # 手动指定顺序：1 → 3 → 2 → 0
)
MarkersprimedACL <- FindMarkers(CER1, group.by = "group", ident.1 = c("CER1_primed_ACL_120h"), ident.2 = c("cRCER1_ACL_CER1+","CER1_primed_DE_72h","cRCER1_DE_72h"),only.pos = T,logfc.threshold = 0.5,min.pct = 0.5)
MarkersnaiveACL <- FindMarkers(ss22411, group.by = "group", ident.1 = c("cRCER1_ACL_CER1+"), ident.2 = c("CER1_primed_ACL_120h"),logfc.threshold = 0.5,min.pct = 0.5)
MarkersprimedDE <- FindMarkers(ss22411, group.by = "group", ident.1 = c("CER1_primed_DE_72h"), ident.2 = c("cRCER1_ACL_CER1+","CER1_primed_ACL_120h","cRCER1_DE_72h"),only.pos = T,logfc.threshold = 0.5,min.pct = 0.5)
MarkersnaiveDE <- FindMarkers(ss22411, group.by = "group", ident.1 = c("cRCER1_DE_72h"), ident.2 = c("CER1_primed_DE_72h","cRCER1_ACL_CER1+","CER1_primed_ACL_120h"),only.pos = T,logfc.threshold = 0.5,min.pct = 0.5)

MarkersnaiveACLvsprimed <- FindMarkers(CER1, group.by = "group", ident.1 = c("cRCER1_ACL_CER1+"), ident.2 = c("CER1_primed_ACL_120h"),logfc.threshold = 0.5,min.pct = 0.5,,pseudocount.use = 0)

# 获取原始 count matrix（不是 normalized）
counts <- GetAssayData(CER1, slot = "counts", assay = "RNA")
counts_cpm <- t(t(counts) / Matrix::colSums(counts)) * 1e6



# 获取 clusterA 的细胞名称
cells_naiveACL <- WhichCells(CER1, idents = "cRCER1_ACL_CER1+")
cells_primedACL <- WhichCells(CER1, idents = "CER1_primed_ACL_120h")
counts_naiveACL <- counts[, cells_naiveACL]
counts_primedACL <- counts[, cells_primedACL]

expr_naiveACL <- Matrix::rowSums(counts_naiveACL > 10)
expr_primedACL <- Matrix::rowSums(counts_primedACL > 10)

# 只要任一 cluster 中满足 ≥24 个细胞 count > 10
expr_filterACL <- (expr_naiveACL >= 24) | (expr_primedACL >= 24)
genes_pass_filter <- names(expr_filterACL[expr_filterACL])

FMarkersnaiveACLvsprimed <- MarkersnaiveACLvsprimed[rownames(MarkersnaiveACLvsprimed) %in% genes_pass_filter, ]
write.table(FMarkersnaiveACLvsprimed,"MarkersnaiveACLvsprimedcount>10.txt",sep = "\t")

MarkersnaiveDEvsprimed <- FindMarkers(CER1, group.by = "group", ident.1 = c("cRCER1_DE_72h"), ident.2 = c("CER1_primed_DE_72h"),logfc.threshold = 0.5,min.pct = 0.5,pseudocount.use = 0)

cells_naiveDE <- WhichCells(CER1, idents = "cRCER1_DE_72h")
cells_primedDE <- WhichCells(CER1, idents = "CER1_primed_DE_72h")
counts_naiveDE <- counts[, cells_naiveDE]
counts_primedDE <- counts[, cells_primedDE]

expr_naiveDE <- Matrix::rowSums(counts_naiveDE > 10)
expr_primedDE <- Matrix::rowSums(counts_primedDE > 10)

# 只要任一 cluster 中满足 ≥24 个细胞 count > 10
expr_filterDE <- (expr_naiveDE >= 24) | (expr_primedDE >= 24)
genes_pass_filterDE <- names(expr_filterDE[expr_filterDE])

FMarkersnaiveDEvsprimed <- MarkersnaiveDEvsprimed[rownames(MarkersnaiveDEvsprimed) %in% genes_pass_filterDE, ]
write.table(FMarkersnaiveDEvsprimed,"MarkersnaiveDEvsprimedcount>10.txt",sep = "\t")


VlnPlot(CER1, features = "nCount_RNA")
VlnPlot(CER1, features = c("TTR","ALDH1A1","APOA4","RBP4","SERPINA1","AFP"),slot = "data")
?FindAllMarkers
DefaultAssay(CER1) <- "RNA"
markerCER1 <- FindAllMarkers(CER1, only.pos = T, logfc.threshold = 0.5,pseudocount.use = 1e-5,min.pct = 0.5)
markerCER1$ratio <- markerCER1$pct.1/markerCER1$pct.2
filtered_markersCER1 <- subset(markerCER1, p_val_adj < 0.05)
write.table(filtered_markers, "250601CERendomarkersFC0.5.txt", sep = "\t")

VlnPlot(ss22411, features = c("HAND1"))
marker.list <- filtered_markers %>% 
  mutate(gene = rownames(.)) %>% 
  group_by(cluster) %>% 
  slice_max(ratio, n = 15) %>% 
  ungroup()
ss22504wt <- NormalizeData(ss22504wt)
"%nin%" = Negate("%in%")
ss22504wt <- subset(ss22504,subset= group %in% c("WT_ACL_D4","WT_DE_D4","cRWT_ACL_D6"))
ss22411noys <- subset(ss22411,subset = group %nin% c("cRCER1_ACL_CER1+_AB","cRCER1_ACL_CER1-_AB","cRH9_ACL_PD+_AB","cRH9_ACL_PD-_AB"))
  VlnPlot(ss22411noys, features = c("AFP",""),group.by = "group",slot = "data")
VlnPlot(ss22411, features = c("AFP","ALDH1A1","VTN","DCN","CER1"),group.by = "group",slot = "data")
VlnPlot(ss22411, features = c("FN1","COL4A1","COL4A2"),group.by = "group",slot = "data")

#########
pdf("250601CERpcaheatmap.pdf", width =6, height = 4)
PCHeatmap(CER1,dims = 1:6,balanced = T,cells = 100)
dev.off()
##
library(destiny)
mat <- as.matrix(CER1@assays$RNA@data)
dm <- DiffusionMap(t(mat))

# 可视化：每个点颜色是 group / cell type
plot(dm$DC1, dm$DC2, col = CER1$group, pch = 16,label = T)
legend("topright",                             # 图例位置
       legend = unique(CER1$group),           # 图例标签
       col = unique(CER1$group),              # 图例颜色（需与上面一致）
       pch = 16,                              # 图例图形
       title = "Group")

DoHeatmap(CER1,features = filtered_markers$gene,slot = "scale.data",size = 4,label = F)+ scale_fill_gradientn(colors = c("#0E53AB", "#eeeeee", "#ef4242"))
DoHeatmap(CER1,features = markersum$gene,slot = "scale.data",size = 4,label = F)+ scale_fill_gradientn(colors = c("#1ECBE1", "#eeeeee", "#E1341E"))
?DoHeatmap
CER1 <- ScaleData(CER1)

DimPlot(CER1,group.by = "group",reduction = "pca")

hypoendo$hypo.ident <- Idents(hypolin)
PCHeatmap(CER1,dims = 1:6,balanced = T,cells = 100)


sc.list <- c(reference, ss22411)

# select features that are repeatedly variable across datasets for integration
?SelectIntegrationFeatures
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 20000)
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, anchor.features = features)
sc.combined <- IntegrateData(anchorset = sc.anchors, k.weight = 30) 
sc.combined <- FindVariableFeatures(sc.combined, selection.method = "vst", nfeatures = 10000)
DefaultAssay(sc.combined) <- "integrated"

sc.combined <- NormalizeData(sc.combined, normalization.method = "RC", scale.factor = 10000)
sc.combined <- ScaleData(sc.combined, verbose = FALSE)
sc.combined <- RunPCA(sc.combined, npcs = 30, verbose = FALSE)
sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindNeighbors(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindClusters(sc.combined, resolution = 0.4)
saveRDS(sc.combined, "ss22411wiCER1endodermwithembryo.RDS")
sc.combined <- readRDS("ss22411wiCER1endodermwithembryo.RDS")


DimPlot(sc.combined,reduction = "umap",group.by = "group")

reference <- readRDS("/Users/hanayoo/Desktop/cam/humansinglecell/CPM7humanplusmamembryoonly240828.RDS")

cells_CER1 <- colnames(CER1)
cells_hypolin <- colnames(hypolin)

# 合并并去重（并集）
cells_union <- union(cells_CER1, cells_hypolin)

# 从 sc.combined 中提取这些细胞
hypoendo <- subset(sc.combined, cells = cells_union)

hypoendo <- ScaleData(hypoendo, verbose = FALSE)
hypoendo <- RunPCA(hypoendo, npcs = 30, verbose = FALSE)
hypoendo <- RunUMAP(hypoendo, reduction = "pca", dims = 1:15)
hypoendo <- FindNeighbors(hypoendo, reduction = "pca", dims = 1:15)
hypoendo <- FindClusters(hypoendo, resolution = 0.4)

hypoendo$hypo.ident <- Idents(hypolin)
hypoendo$hypo.ident <- as.character(hypoendo$hypo.ident)
hypoendo$hypo.ident[is.na(hypoendo$hypo.ident)] <- hypoendo$group[is.na(hypoendo$hypo.ident)]
DimPlot(hypoendo,reduction = "pca",group.by = "hypo.ident")



MarkersnaiveACLvsprimed <- FindMarkers(CER1, group.by = "group", ident.1 = c("cRCER1_ACL_CER1+"), ident.2 = c("CER1_primed_ACL_120h"),pseudocount.use = 1e-5)
MarkersnaiveDEvsprimed <- FindMarkers(CER1, group.by = "group", ident.1 = c("cRCER1_DE_72h"), ident.2 = c("CER1_primed_DE_72h"),pseudocount.use = 0)

write.table(MarkersnaiveACLvsprimed,"MarkersnaiveACLvsprimedtotal.txt",sep = "\t")
write.table(MarkersnaiveDEvsprimed,"MarkersnaiveDEvsprimedtotal.txt",sep = "\t")
