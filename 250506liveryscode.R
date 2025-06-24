library(dplyr)
library(pheatmap)


SYSinvitro <- subset(ss22411, subset = group %in% c("cRCER1_ACL_CER1+_AB","cRH9_ACL_PD+_AB"))
SYSinvitro$time2 <- SYSinvitro$group

SYSearly <- readRDS("/Users/hanayoo/Desktop/cam/YSseq/E-MTAB-10888/YSss2count/YSearly.RDS")
SYSearly <- subset(SYSearly, subset = species == "human")
SYSearly$time2 <- "YSE_<3W"
DefaultAssay(SYSearly) <- "RNA"


SYSlate <- readRDS("/Users/hanayoo/Desktop/cam/YSseq/E-MTAB-10888/YSss2count/YSlate.RDS")
SYSlate$time2<-SYSlate$time
metadata <- read.table("/Users/hanayoo/Desktop/cam/YSseq/E-MTAB-10888/YSss2count/metadata.txt", header = T)
SYSlate <- AddMetaData(SYSlate, metadata = metadata, col.name = c("cell.labels","cluster","stage","week","Day2"))
SYSlate$time2 <- paste0("YSE_", SYSlate$Day2)

liver <- readRDS("/Users/hanayoo/Desktop/cam/liverseq/liver-V4.RDS")
special_rename <- c(
  "COX1" = "MT-CO1",
  "COX2" = "MT-CO2",
  "COX3" = "MT-CO3",
  "CYTB" = "MT-CYB"
)

# 其他线粒体基因需要加前缀“MT-”，如 ND1 -> MT-ND1
# 这些是标准线粒体基因名（未带MT-前缀的）
basic_mito_genes <- c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8", "CO1", "CO2", "CO3", "CYB")

# 提取 counts 和 data 矩阵
counts <- liver[["RNA"]]@counts
data <- liver[["RNA"]]@data

# 第一步：处理特殊名称合并
for (old in names(special_rename)) {
  new <- special_rename[[old]]
  
  if (old %in% rownames(counts)) {
    if (new %in% rownames(counts)) {
      counts[new, ] <- counts[new, ] + counts[old, ]
      data[new, ] <- data[new, ] + data[old, ]
    } else {
      rownames(counts)[rownames(counts) == old] <- new
      rownames(data)[rownames(data) == old] <- new
    }
  }
}

# 第二步：处理无MT-前缀的基因（加MT-并合并）
for (base in basic_mito_genes) {
  old <- base
  new <- paste0("MT-", base)
  
  if (old %in% rownames(counts)) {
    if (new %in% rownames(counts)) {
      counts[new, ] <- counts[new, ] + counts[old, ]
      data[new, ] <- data[new, ] + data[old, ]
    } else {
      rownames(counts)[rownames(counts) == old] <- new
      rownames(data)[rownames(data) == old] <- new
    }
  }
}

# 最后：确保无重复并清理原始行名
counts <- counts[!duplicated(rownames(counts)), ]
data <- data[!duplicated(rownames(data)), ]

# 更新 Seurat 对象
liver[["RNA"]]@counts <- counts
liver[["RNA"]]@data <- data


sc.list <- c(liver,SYSearly,SYSlate,SYSinvitro)
sc.list <- lapply(X = sc.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "RC", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 50000)
})
p2 <-DimPlot(liver)
p3 <-DimPlot(liver, group.by = "SAMPLE")
p1 <- FeaturePlot(liver, features = c("AFP","ALB","TTR"),max.cutoff = 50)
p1+p2+p3

# select features that are repeatedly variable across datasets for integration
?SelectIntegrationFeatures
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 20000)
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, anchor.features = features)
sc.combined <- IntegrateData(anchorset = sc.anchors, k.weight = 30) 
sc.combined <- FindVariableFeatures(sc.combined, selection.method = "vst", nfeatures = 50000)
DefaultAssay(sc.combined) <- "integrated"

meta <- read.table("totalmeta2.txt", sep = "\t", header = TRUE, row.names = 1)
endoderm <- AddMetaData(object = endoderm, metadata = meta, col.name = c("time","time2","type"))
sc.combined <- AddMetaData(object = sc.combined, metadata = meta, col.name = c("time","time2","type"))

sc.combined <- NormalizeData(sc.combined, normalization.method = "RC", scale.factor = 10000)
sc.combined <- ScaleData(sc.combined, verbose = FALSE)
sc.combined <- RunPCA(sc.combined, npcs = 30, verbose = FALSE)
sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindNeighbors(sc.combined, reduction = "pca", dims = 1:15)
sc.combined <- FindClusters(sc.combined, resolution = 0.4)

meta <- read.table("/Users/hanayoo/Desktop/cam/YSLIVER/totalmeta2.txt", sep = "\t", header = TRUE, row.names = 1)
endoderm <- AddMetaData(object = endoderm, metadata = meta, col.name = c("time","time2","type"))
sc.combined <- AddMetaData(object = sc.combined, metadata = meta, col.name = c("time","time2","type"))

saveRDS(sc.combined, "250505YSliverinvitro.RDS")
sc.combined <- readRDS("/Users/hanayoo/Desktop/cam/seq2504-mt/250505YSliverinvitro.RDS")

P1 <- DimPlot(sc.combined,group.by = "time2",reduction = "pca")
P2 <-DimPlot(sc.combined,group.by = "cluster",reduction = "pca")
P3 <-FeaturePlot(sc.combined,features = c("AFP"),slot = "data",min.cutoff = 0,,max.cutoff = 50,reduction = "pca")
P3 <- FeaturePlot(sc.combined,features = c("TTR"),slot = "data",min.cutoff = 0,max.cutoff = 50,reduction = "pca")
P4 <- DimPlot(sc.combined,label = T,label.size = 2,reduction = "pca")
p <- wrap_plots(P1,P2,P3,P4,ncol = 2)
p

endoderm <- subset(sc.combined, idents = c(0,4,7,8,9))
endoderm <- NormalizeData(endoderm, normalization.method = "RC", scale.factor = 10000)
endoderm <- ScaleData(endoderm, verbose = FALSE)
endoderm <- RunPCA(endoderm, npcs = 30, verbose = FALSE)
endoderm <- RunUMAP(endoderm, reduction = "pca", dims = 1:15)
endoderm <- FindNeighbors(endoderm, reduction = "pca", dims = 1:15)
endoderm <- FindClusters(endoderm, resolution = 0.4)

DimPlot(endoderm,group.by = "time2",reduction = "pca")
pdf("250505_ysliverendoderm.pdf", width =6, height = 4)
DimPlot(endoderm,group.by = "time2",reduction = "pca")
dev.off()

pdf("250505_ysliversombined.pdf", width =12, height = 8)
plot(p)
dev.off()


pdf("250505_yslivermarkers.pdf", width =12, height = 4)
VlnPlot(endoderm, features = c("AFP","ALB","FTH1","CYP3A4"),group.by = "time2",slot = "data")
dev.off()
VlnPlot(endoderm, features = c("APOA1-AS","CFH","ABCB9","HAMP","MYL5","CFHR1","NR1H4","LIN28B"),group.by = "time2",slot = "data", ncol = 4)
DefaultAssay(endoderm) <- "RNA"

Idents(endoderm) <- endoderm$time2
saveRDS(endoderm, "LIVERYSendoderm.RDS")
endoderm <- readRDS("/Users/hanayoo/Desktop/cam/seq2504-mt/LIVERYSendoderm.RDS")
pca.stdev <- endoderm@reductions$pca@stdev  # 各PC的标准差
pca.var <- pca.stdev^2                 # 方差 = stdev^2
var.percent <- pca.var / sum(pca.var) * 100 
var.percent


invivoendo <- subset(endoderm, subset = time2 %in% c("liver_ADULT","liver_FETAL","YSE_<3W","YSE_5W","YSE_7W"))
markers <- FindAllMarkers(invivoendo,logfc.threshold = 1, group.by = "time2",only.pos = T,min.pct = 0.5)
markers <- markers %>%
  mutate(ratio = pct.1 / pct.2)

write.table(markers, "allmarkers.txt",sep = "\t")

gene.lists <- lapply(sc.list, function(obj) rownames(obj[["RNA"]]@counts))

# 2. 取交集，得到四个对象共有的基因
common.genes <- Reduce(intersect, gene.lists)

filtered_markers <- markers %>%
  filter(gene %in% common.genes)

top10 <- filtered_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = ratio, n = 15, with_ties = FALSE)

liver10 <- filtered_livermarkers %>%
  slice_max(order_by = ratio, n = 15, with_ties = FALSE)

# 合并两个结果
top_bottom_markers <- bind_rows(top10, liver10)
#endoderm <- ScaleData(endoderm)
marker_genes <- top_bottom_markers %>% pull(gene) %>% unique()

marker_genessort <- c(
  "MFSD2A", "SLPI", "NNMT", "CES1", "CXCL2", "HSD17B13", "SDS", "SLC22A1",
  "CYP2C9", "GRAMD4","C4BPA","MYOM1", "APOA1-AS", "CFH",  "HAMP", "MYL5", 
  "ABCB9", "CFHR1", "NR1H4", "COL2A1", "NDN", "VCAM1", "CYP19A1", "LDOC1", "DLK1",
  "STARD4-AS1", "BEX2", "CD4", "RBP1", "LIN28B", "DDX47", "LAMA1", "NDUFA13",
  "CTSV", "MRPL38", "CLDN6", "ANGPT1", "WNT11", "RPS10", 
  "LIPH", "MSR1", "ALOX15B", "PRDM1", "RGS16", "TMEM144", "FKBP10", "SLC5A11",
  "GNG11","MCOLN3", "AXL",
  "COL3A1", "PARM1", "PDGFA", "SLIT2", "EFNA5", "CBX6", "F3", "SLIT3"
)

marker_genessort <- c(
  "MFSD2A", "SLPI", "NNMT", "CES1", "CXCL2", "APOA1-AS", "CFH", "MYL5", 
  "ABCB9", "COL2A1", "NDN", "VCAM1", "CYP19A1", "LDOC1", "DLK1",
  "LIN28B", "DDX47", "LAMA1", "NDUFA13",
  "CTSV", "MRPL38", "CLDN6", "ANGPT1", "RPS10", 
  "LIPH", "MSR1", "ALOX15B", "PRDM1", "RGS16", "TMEM144", "FKBP10", "SLC5A11",
  "MCOLN3", "AXL",
  "COL3A1", "PARM1", "PDGFA", "SLIT2", "EFNA5", "CBX6"
)


DefaultAssay(endoderm) <- "RNA"
VlnPlot(endoderm, features = c("LIN28B", "DDX47", "LAMA1", "NDUFA13",
                               "CTSV", "MRPL38", "CLDN6", "ANGPT1"),ncol = 4)
DoHeatmap(endoderm, features = marker_genessort, size = 3,label = F) +  
  scale_fill_gradientn(colors = c("#0E53AB", "#FFFFFF", "#EF4242"))
pdf("250619_yslivermarkersheatmapresort.pdf", width =5, height = 6)
DoHeatmap(endoderm, features = marker_genessort, size = 3,label = F) +  
 scale_fill_gradientn(colors = c("#0E53AB", "#FFFFFF", "#EF4242"))
dev.off()




markersliver <- FindMarkers(invivoendo,logfc.threshold = 1,only.pos = T,min.pct = 0.5,ident.1 = c("liver_ADULT","liver_FETAL"),ident.2 = c("YSE_<3W","YSE_5W","YSE_7W"))
markersliver <- markersliver %>%
  mutate(ratio = pct.1 / pct.2)
markersliver$gene <- row.names(markersliver)
filtered_livermarkers <- markersliver %>%
  filter(gene %in% common.genes)


write.table(markersliver, "markersliver.txt",sep = "\t")

# 提取前3个主成分的坐标
pca_embed <- Embeddings(endoderm, reduction = "pca")[, 1:3]

# 创建 dataframe，并加入 metadata 中你关心的列（例如 mesosub$group）
plot_df <- data.frame(pca_embed, 
                      group = endoderm$time2,  # 可换成你想分组的列，比如 mesosub$seurat_clusters
                      cell = colnames(endoderm))

plot_ly(plot_df,
        x = ~PC_1, y = ~PC_2, z = ~PC_3,
        color = ~group,
        text = ~cell,
        type = "scatter3d",
        mode = "markers" ,
        marker = list(size = 5)) %>%  # 👈 调小点的大小（默认大约是 5）%>%
  layout(scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")
  ))

Idents(endoderm) <- "time2"

counts <- GetAssayData(endoderm, slot = "counts", assay = "RNA")

# 取每个 group 的 pseudobulk
group_ids <- Idents(endoderm)
group_levels <- levels(group_ids)

# 只保留感兴趣的基因
common_genes <- intersect(rownames(counts), marker_genessort)
counts <- counts[common_genes, ]

# 创建 pseudobulk 表达矩阵（每列一个 group）
pseudobulk_mat <- sapply(group_levels, function(g) {
  cells_in_group <- WhichCells(endoderm, idents = g)
  if (length(cells_in_group) == 0) return(rep(0, length(common_genes)))
  Matrix::rowSums(counts[, cells_in_group, drop = FALSE])
})

rownames(pseudobulk_mat) <- common_genes

# log1p 转换
log_mat <- log1p(pseudobulk_mat)

# row scale（按基因标准化）
scaled_mat <- t(scale(t(log_mat)))
scaled_mat <- scaled_mat[marker_genessort, , drop = FALSE]
pheatmap(scaled_mat,
         cluster_rows = F,
         cluster_cols = F,
         fontsize_row = 8,
         fontsize_col = 10,
         color = colorRampPalette(c("#0E53AB", "#FFFFFF", "#EF4242"))(100))
saveRDS(common.genes,"liveryscommongenes.RDS")

endoderm$time2 <- factor(endoderm$time2, levels = rev(c("cRCER1_ACL_CER1+_AB", "cRH9_ACL_PD+_AB","liver_ADULT","liver_FETAL","YSE_<3W","YSE_5W", 
                                                      "YSE_7W")))
Idents(endoderm) <- endoderm$time2
# 使用 DotPlot 直接画图
dot_data <- DotPlot(endoderm, features = marker_genessort, group.by = "time2")$data

dot_data <- dot_data %>%
  mutate(log_avg_exp = log1p(avg.exp)) %>%
  group_by(features.plot) %>%
  mutate(scaled_log_avg_exp = rescale(log_avg_exp))

ggplot(dot_data, aes(x = features.plot, y = id, size = pct.exp, color = scaled_log_avg_exp)) +
  geom_point() +
  scale_color_gradient(low = "grey", high = "#0e53ab") +
  theme_minimal() +
  labs(x = "Gene", y = "Group", size = "% Expressing", color = "Scaled Avg Exp") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(scuttle)
library(pheatmap) 
pseudobulk_expression_matrix <- Seurat::AverageExpression(endoderm,
                                                          group.by = "time2", # Use the name of your metadata column
                                                          assays = "RNA",
                                                          slot = "counts" # CRITICAL: Get raw counts
)$RNA # Extract the matrix for the "RNA" assay

endoderm <- FindVariableFeatures(endoderm, selection.method = "vst", nfeatures = 2000)

# 获取 HVGs 列表
hvgs_from_endoderm <- VariableFeatures(endoderm)
final_gene_list <- intersect(common.genes, hvgs_from_endoderm)
filtered_pseudobulk_counts <- pseudobulk_expression_matrix[final_gene_list, ]
total_counts_per_sample <- colSums(filtered_pseudobulk_counts)
pseudobulk_cpm <- sweep(filtered_pseudobulk_counts,
                        MARGIN = 2,
                        STATS = total_counts_per_sample / 1000000,
                        FUN = "/")
correlation_matrix_cpm <- cor(pseudobulk_cpm, method = "spearman") # Or "spearman"
my_colors_1 <- colorRampPalette(c("#0e53ab", "#ffffff", "#ef4242"))(100)
pheatmap(correlation_matrix_cpm,
         display_numbers = TRUE,
         cluster_rows = FALSE, # 关闭行聚类
         cluster_cols = FALSE, # 关闭列聚类
         main = "Gene Expression Correlation (common genes, CPM)")

pdf("250610_correlationliverendodermspearmanhvgs200022.pdf", width =5, height = 5)
pheatmap(correlation_matrix_cpm,
         display_numbers = TRUE,
         cluster_rows = FALSE, # 关闭行聚类
         cluster_cols = FALSE, # 关闭列聚类
         color = my_colors_1,
         main = "Gene Expression Correlation (common genes, CPM)")
         
dev.off()



