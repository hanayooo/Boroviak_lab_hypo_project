rm(list = ls())
setwd("/Users/hanayoo/Desktop/cam/seq2504-mt")

library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(ggplot2)

# 读取 count matrix（假设是 tab 分隔）
df <- read.delim("merged_counts.txt", header = TRUE, stringsAsFactors = FALSE,check.names = FALSE)

# 拆分第一列：Geneid.gene_name -> geneid + genename
gene_info <- do.call(rbind, strsplit(df[[1]], "\\*"))
colnames(gene_info) <- c("geneid", "genename")

# 添加新列并删除原来的
df <- cbind(gene_info, df[ , -1])

clean_colnames <- function(x) {
  if (grepl("SLX-25104\\.", x)) {
    # 提取 i7xx-i5yy 片段
    matches <- regmatches(x, regexpr("i7[^.]+-i5[^.]+", x))
    return(matches)
  } else {
    return(x)  # geneid / genename 不动
  }
}

# 更新列名
colnames(df) <- sapply(colnames(df), clean_colnames)

# 保存新表格
write.table(df, "cleaned_count_matrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 删除 gene_name 为 "gene_name" 的行
df <- df[df$genename != "gene_name", ]
df$geneid <- NULL
df[ , -1] <- lapply(df[ , -1], as.numeric)

df$total_count <- rowSums(df[ , -1])

# 对 genename 去重：保留总和最高那一行
df <- df[order(-df$total_count), ]                  # 按 total_count 降序
df <- df[!duplicated(df$genename), ]                # 保留第一次出现的 genename（最大值）
df$total_count <- NULL                              # 删除辅助列

# 设置 genename 为行名
rownames(df) <- df$genename
df$genename <- NULL
write.table(df, "cleaned_count_matrix_rmdup.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 读入 barcode 映射表
barcode_map <- read.table("barcodes.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# 构建映射命名向量
barcode_dict <- setNames(barcode_map$sample, barcode_map$barcodes)


# 替换列名
colnames(df) <- ifelse(colnames(df) %in% names(barcode_dict),
                       barcode_dict[colnames(df)],
                       colnames(df))

# 保存结果
write.csv(df, "2411_count_matrix_renamed_mt.csv")
df <- read.csv("2411_count_matrix_renamed_mt.csv",header = T,row.names = 1)
ss22504 <- CreateSeuratObject( counts = df, min.cells = 1, min.features = 200)
meta <- read.table("metadata.txt", sep = "\t", header = TRUE, row.names = 1)
ss22504 <- AddMetaData(object = ss22504, metadata = meta, col.name = c("barcodes","group") )

"%nin%" = Negate("%in%")
ss22504 <- subset(ss22504, subset = group %nin% c("CM_J") )
###ss22504 <- subset(ss22504, subset = group %nin% c("MS615_3D","MS615_2D") )

ss22504 <- FindVariableFeatures(ss22504, selection.method = "vst", nfeatures = 60000)
ss22504 <- NormalizeData(ss22504, normalization.method = "RC", scale.factor = 10000)
ss22504 <- ScaleData(ss22504, verbose = FALSE)
ss22504 <- RunPCA(ss22504, npcs = 30, verbose = FALSE)
ss22504 <- RunUMAP(ss22504, reduction = "pca", dims = 1:15)
ss22504 <- FindNeighbors(ss22504, reduction = "pca", dims = 1:15)
ss22504 <- FindClusters(ss22504, resolution = 0.4)

MSExMes<- subset(ss22504,subset = group %in% c("MS615_3D","MS615_2D"))
mergedss22411 <- merge(ss22411, y = MSExMes, add.cell.ids = c("A", "B"), project = "mergedss22411")

mergedss22411 <- FindVariableFeatures(mergedss22411, selection.method = "vst", nfeatures = 60000)
mergedss22411 <- NormalizeData(mergedss22411, normalization.method = "RC", scale.factor = 10000)
mergedss22411 <- ScaleData(mergedss22411, verbose = FALSE)
mergedss22411 <- RunPCA(mergedss22411, npcs = 30, verbose = FALSE)
mergedss22411 <- RunUMAP(mergedss22411, reduction = "pca", dims = 1:15)
mergedss22411 <- FindNeighbors(mergedss22411, reduction = "pca", dims = 1:15)
mergedss22411 <- FindClusters(mergedss22411, resolution = 0.4)

saveRDS(mergedss22411,"ss22411v4wiMSExMes")
DimPlot(mergedss22411, group.by = "group")
DimPlot(ss22504, group.by = "group")

sc.list <- c(reference, mergedss22411)

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
saveRDS(sc.combined,"ss22411v4wiMSExMeswithref")
sc.combined$meso.ident <- sc.combined$orianos2
mesoderm <- c("Advanced Mesoderm","Axial Mesoderm","cRH9_ACL_PD-",
              "cRH9_ACL_PD-_AB","Emergent Mesoderm","ExE Mesoderm",
              "MS615_2D","MS615_3D","Nascent Mesoderm","Primitive Streak",
              "cRCER1_ACL_CER1-","cRCER1_ACL_CER1-_AB")
exmes <-c("Advanced Mesoderm","cRH9_ACL_PD-",
          "cRH9_ACL_PD-_AB","ExE Mesoderm",
          "MS615_2D","MS615_3D",
          "cRCER1_ACL_CER1-","cRCER1_ACL_CER1-_AB")

# 对 meso.ident 为 NA 的细胞，用 group 的值填充
sc.combined$meso.ident[is.na(sc.combined$meso.ident)] <- sc.combined$group[is.na(sc.combined$meso.ident)]
DimPlot(sc.combined, group.by = "meso.ident")


mesosub <-subset(sc.combined, subset = meso.ident %in% mesoderm)
mesosub <- ScaleData(mesosub)
mesosub <- RunPCA(mesosub, npcs = 30, verbose = FALSE)
mesosub <- RunUMAP(mesosub, reduction = "pca", dims = 1:15)
mesosub <- FindNeighbors(mesosub, reduction = "pca", dims = 1:15)
mesosub <- FindClusters(mesosub, resolution = 0.4)
DimPlot(mesosub, group.by = "meso.ident",reduction = "pca",dims = c(1,3))

# 提取前3个主成分的坐标
pca_embed <- Embeddings(mesosub, reduction = "pca")[, 1:3]

# 创建 dataframe，并加入 metadata 中你关心的列（例如 mesosub$group）
plot_df <- data.frame(pca_embed, 
                      group = mesosub$meso.ident.sort,  # 可换成你想分组的列，比如 mesosub$seurat_clusters
                      cell = colnames(mesosub))

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
write.table(plot_df,"PCA-3Dmeso.txt",sep = "\t")

VlnPlot(mesosub, features = c("COL3A1","POSTN","LUM","DCN","HAND1","FOXF1"),group.by = "meso.ident.sort",slot = "data",label)
mesosub <- FindVariableFeatures(mesosub, selection.method = "vst", nfeatures = 60000)
DefaultAssay(mesosub) <-"RNA"

pdf("250506_exmesmarkers.pdf", width =12, height = 9)
VlnPlot(mesosub, features = c("COL3A1","POSTN","LUM","DCN","HAND1","FOXF1"),group.by = "meso.ident.sort",slot = "data")
dev.off()


# 复制原始值
mesosub$meso.ident.sort <- mesosub$meso.ident

# 替换指定分类名
mesosub$meso.ident.sort[mesosub$meso.ident == "Axial Mesoderm"] <- "N_Axial Mesoderm"
mesosub$meso.ident.sort[mesosub$meso.ident == "Emergent Mesoderm"] <- "N_Emergent Mesoderm"

saveRDS(mesosub,"mesosub.RDS")


