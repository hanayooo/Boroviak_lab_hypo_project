rm(list = ls())
setwd("/Users/hanayoo/Desktop/cam/seq2411-mt")
#install.packages("Seurat")

library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

ss22411count<- read.table("merged_counts_name.txt",sep = "\t",header = T,nrows = 58884,row.names = 1,check.names = FALSE)

#####
rownames_split <- strsplit(rownames(ss22411count), "\\*")
genenames <- sapply(rownames_split, function(x) x[2]) 
ss22411count$genename <- genenames
ss22411count$row_id <- seq_len(nrow(ss22411count))
ss22411count <- ss22411count %>%
  group_by(genename) %>%
  arrange(desc(H9_primed_1), row_id) %>%  # 按排序列降序，行号升序排列
  slice(1) %>%  # 只保留分组内第一行
  ungroup()
ss22411count <- as.data.frame(ss22411count)
rownames(ss22411count) <- ss22411count$genename

ss22411count$genename <- NULL
ss22411count$row_id <- NULL

write.csv(ss22411count, "250505-ss22411countrmdup-mt.csv")

ss22411count <- read.csv("/Users/hanayoo/Desktop/cam/seq2411-mt/250505-ss22411countrmdup-mt.csv",row.names = 1,check.names = F)
######

ss22411 <- CreateSeuratObject(counts = ss22411count, project = "ss22411", min.cells = 3, min.features = 1000)

ss22411 <- FindVariableFeatures(ss22411, selection.method = "vst", nfeatures = 60000)

ss22411 <- NormalizeData(ss22411, normalization.method = "RC", scale.factor = 10000)
ss22411 <- ScaleData(ss22411, verbose = FALSE)
ss22411 <- RunPCA(ss22411, npcs = 30, verbose = FALSE)
ss22411 <- RunUMAP(ss22411, reduction = "pca", dims = 1:15)
ss22411 <- FindNeighbors(ss22411, reduction = "pca", dims = 1:15)
ss22411 <- FindClusters(ss22411, resolution = 1)

####
sample <- colnames(ss22411count) 
groups <- sub("_\\d+$", "", sample)
df <- data.frame(cells = sample, group = groups, stringsAsFactors = FALSE)
####

ss22411 <- AddMetaData(object = ss22411, metadata = df, col.name = c("group"))

p1 <- DimPlot(ss22411, reduction = "umap", group.by = "group",label = T, label.size = 2, pt.size = 0.2)
p2 <- DimPlot(ss22411, reduction = "umap",label = T, label.size = 2, pt.size = 0.2)
DimPlot(ss22411, reduction = "pca",label = T, label.size = 2, pt.size = 0.2)
FeaturePlot(ss22411, features = c("TTR","AFP","FABP1","ALB"))
p1+p2

saveRDS(ss22411, "ss22411-1123-v4-mt.RDS")
ss22411 <- readRDS("ss22411-1123-v4-mt.RDS")

ss22411$group <- factor(ss22411$group, levels = c("CER1_PXGL","H9_PXGL","CER1_primed","H9_primed", 
                                                      "cRCER1_ACL_CER1+", "cRH9_ACL_PD+","cRCER1_ACL_CER1-", "cRH9_ACL_PD-", 
                                                      "cRCER1_ACL_CER1+_AB", "cRH9_ACL_PD+_AB","cRCER1_ACL_CER1-_AB", "cRH9_ACL_PD-_AB"))
Idents(ss22411) <- ss22411$group
# 使用 DotPlot 直接画图
VlnPlot(ss22411, features = c("POU5F1","KLF4","SOX17","LEFTY1","AFP","TTR","HAND1","FOXF1"),ncol = 4)



######remove cells for DE/VE comparison####
"%nin%" = Negate("%in%")
ss22411 <- subset(ss22411, subset = group %nin% c("CER1_primed_DE_72h","cRCER1_DE_72h","CER1_primed_ACL_120h","CER1_primed_ACL_72h"))
DimPlot(ss22411,group.by = "group",label = T)

reference <- readRDS("/Users/hanayoo/Desktop/cam/humansinglecell/CPM7humanplusmamembryoonly240828.RDS")
reference$rough.ident <- Idents(reference)

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

p1 <- DimPlot(sc.combined, reduction = "umap",group.by = "group",label = F,label.size = 3)
#p2 <- DimPlot(sc.combined, reduction = "umap",group.by = "active.ident",label = T,label.size = 4) + NoLegend()+coord_cartesian(xlim = c(5，5), ylim = c(0, 10))
p3<- DimPlot(sc.combined, reduction = "umap",group.by = "orianos2",label = F,label.size = 3)
p1+p3

sc.combined$rough.ident[is.na(sc.combined$rough.ident)] <- sc.combined$group[is.na(sc.combined$rough.ident)]

DimPlot(sc.combined, reduction = "umap",group.by = "rough.ident",label = F,label.size = 3)+DimPlot(sc.combined, reduction = "umap",group.by = "orianos2",label = F,label.size = 3)
######
pdf("250606_ss22411vlnplots.pdf", width =14, height = 8)
VlnPlot(ss22411, features = c("POU5F1","KLF4","SOX17","LEFTY1","AFP","TTR","HAND1","FOXF1"),ncol = 4)
dev.off()
#####
hypolin <- readRDS("/Users/hanayoo/Desktop/cam/humansinglecell/hypolin0828.RDS")
meta <- read.table("/Users/hanayoo/Desktop/cam/humansinglecell/cells.txt", sep = "\t", header = TRUE, row.names = 1)
hypolin <- AddMetaData(object = hypolin, metadata = meta, col.name = c("author","oriano","orianos","day","seq","orianos2","subcluster","species","day2","annos") )

hypolin <- FindClusters(hypolin, resolution = 2)
new.cluster.ids<-c("Endo", "Pre_HYPO", "Endo", "Pre_HYPO", "SYSE","Pre_HYPO","VE","Endo","VE","Endo","SYSE"
)
names(new.cluster.ids)<-levels(hypolin)
hypolin<-RenameIdents(hypolin, new.cluster.ids)
DimPlot(hypolin, reduction = "pca",pt.size = 0.2,label = F)


hypo_ids <- Idents(hypolin)
sc.combined$hypo.ident <- NA
sc.combined$hypo.ident[names(hypo_ids)] <- as.character(hypo_ids)
na_cells <- is.na(sc.combined$hypo.ident)
sc.combined$hypo.ident[na_cells] <- sc.combined$group[na_cells]


#####
pdf("250603ss22411dotplot.pdf", width =8, height = 3.5)
ggplot(dot_data, aes(x = features.plot, y = id, size = pct.exp, color = scaled_log_avg_exp)) +
  geom_point() +
  scale_color_gradient(low = "grey", high = "#0e53ab") +
  theme_minimal() + # Start with minimal theme
  labs(x = "Gene", y = "Group", size = "% Expressing", color = "Scaled Avg Exp") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank() , # Remove minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) # 添加黑色边框
  )
dev.off()
####
sc.combined$exmes <- NA  # 先创建列并赋值为NA

# 条件1：orianos2 在指定类别中的，exmes 取 orianos2 的值
sc.combined$exmes[sc.combined$orianos2 %in% c("ExMes_CS5", "ExMes_CS6", "ExMes_CS7", "ExE Mesoderm", "Advanced Mesoderm")] <-
  sc.combined$orianos2[sc.combined$orianos2 %in% c("ExMes_CS5", "ExMes_CS6", "ExMes_CS7", "ExE Mesoderm", "Advanced Mesoderm")]

# 条件2：orianos2 为 NA 的，exmes 取 group 的值
sc.combined$exmes[is.na(sc.combined$orianos2)] <- sc.combined$group[is.na(sc.combined$orianos2)]
p1 <- DimPlot(sc.combined, reduction = "umap",group.by = "exmes",label = F,label.size = 3) + 
  xlim(-4, 4) + 
  ylim(3, 11)
p2 <-DimPlot(sc.combined, reduction = "umap",group.by = "orianos2",label = F,label.size = 3) + 
  xlim(-4, 4) + 
  ylim(3, 11)
p1+p2
pdf("250506_exmesidentorianos2.pdf", width =16, height = 4)
plot(p1+p2)
dev.off()

saveRDS(sc.combined,"seq2411withembryo.RDS")
sc.combined <- readRDS("seq2411withembryo.RDS")

DimPlot(endoderm,group.by = "author")
primed <- subset(sc.combined, subset = group %in% c("CER1_primed","H9_primed"))
write.table(primed@assays[["RNA"]]@counts,"Primedcounts250506.txt",sep = "\t")

sc.combined@meta.data$total_ident <- ifelse(
  is.na(sc.combined@meta.data$orianos2), # 条件：orianos2 列的值是否为 NA
  sc.combined@meta.data$group,         # 如果条件为 TRUE (是 NA)，则使用 group 列的值
  sc.combined@meta.data$orianos2       # 如果条件为 FALSE (不是 NA)，则使用 orianos2 列的值
)

DimPlot(sc.combined, group.by = "total_ident")
"%nin%" <-Negate("%in%")
epiblast <- FindVariableFeatures(epiblast)
epiblast <- subset(sc.combined, subset = total_ident %in% c("CER1_primed","CER1_PXGL","EPI","H9_primed","H9_PXGL","ICM"))
epiblast$author[is.na(epiblast$author)] <- "this study"
epiblast <- subset(epiblast, subset = author %nin% c("Zhou"))
epiblast <- ScaleData(epiblast)
epiblast <- RunPCA(epiblast, npcs = 30, verbose = FALSE)
ElbowPlot(epiblast)
epiblast <- RunUMAP(epiblast, reduction = "pca", dims = 1:10)
epiblast <- FindNeighbors(epiblast, reduction = "pca", dims = 1:15)
epiblast <- FindClusters(epiblast, resolution = 0.4)

epiblast$sorted.date <- NA  # 初始化新 metadata 列
epiblast$sorted.date[epiblast$day2 %in% c("D05", "D05.5", "D06", "D06.5")] <- "D5-D6.5"
epiblast$sorted.date[epiblast$day2 %in% c("D07", "D07.5", "D08")] <- "D7-D8"
epiblast$sorted.date[epiblast$day2 %in% c("D09", "D10", "D11")] <- "D9-D11"
epiblast$sorted.date[epiblast$day2 %in% c("D12", "D13", "D14")] <- "E12-D14"
epiblast$sorted.date[epiblast$day2 %in% c("D16-19")] <- "E16-19"
epiblast$sorted.date[is.na(epiblast$sorted.date)] <- epiblast$group[is.na(epiblast$sorted.date)] 
epiblast <- AddMetaData(object = epiblast, metadata = meta, col.name = c("author","oriano","orianos","day","seq","orianos2","subcluster","species","day2","annos") )
FeaturePlot(epiblast, features = c("KLF4"),min.cutoff = 0, max.cutoff = 4)
DimPlot(epiblast,group.by  = "total_ident",reduction = "umap") 
DimPlot(epiblast,group.by  = "day2",reduction = "umap")
DimPlot(epiblast,group.by  = "author",reduction = "umap") 
VlnPlot(epiblast, group.by = "sorted.date",features = c("SUSD2","OTX2","SOX11","LEFTY1","LIN28A","TFAP2C"))
VlnPlot(epiblast, group.by = "sorted.date",features = c("SUSD2"))
DefaultAssay(epiblast) <- "RNA"

pca_embed <- Embeddings(ss22411, reduction = "pca")[, 1:3]

# 创建 dataframe，并加入 metadata 中你关心的列（例如 mesosub$group）
plot_df <- data.frame(pca_embed, 
                      group = ss22411$seurat_clusters,  # 可换成你想分组的列，比如 mesosub$seurat_clusters
                      cell = colnames(ss22411))
library(plotly)
p <- plot_ly(plot_df,
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
for(i in 1:nrow(plot_df)) {
  p <- p %>% add_trace(
    x = c(plot_df$PC_1[i], plot_df$PC_1[i]),
    y = c(plot_df$PC_2[i], plot_df$PC_2[i]),
    z = c(0, plot_df$PC_3[i]),  # 从Z=0(XY平面)到数据点
    type = "scatter3d",
    mode = "lines",
    line = list(color = "gray", width = 1, dash = "dot"),
    showlegend = FALSE,
    hoverinfo = "none",
    inherit = FALSE  # 不继承颜色等属性
  )
} 
p <- p %>% add_trace(
  x = ~PC_1,
  y = ~PC_2,
  z = rep(0, nrow(plot_df)),  # 所有点在Z=0平面
  type = "scatter3d",
  mode = "markers",
  color = ~group,
  marker = list(size = 3, opacity = 0.7),
  text = ~cell,
  name = "XY投影",
  inherit = FALSE
)


write.table(plot_df,"250610-PCA-3Dss22411.txt",sep = "\t")


# 假设你已经有一个 Seurat 对象，例如：seu

library(scales)
# 指定要可视化的一组基因
genes.to.plot <- c("POU5F1", "NANOG", "KLF4", "KLF17", "GATA6","SOX17","LEFTY1","TTR","APOB","AFP","HAND1","FOXF1","POSTN","DCN")
ss22411$group <- factor(ss22411$group, levels = rev(c("CER1_PXGL","H9_PXGL","CER1_primed","H9_primed", 
                                                  "cRCER1_ACL_CER1+", "cRH9_ACL_PD+","cRCER1_ACL_CER1-", "cRH9_ACL_PD-", 
                                                  "cRCER1_ACL_CER1+_AB", "cRH9_ACL_PD+_AB","cRCER1_ACL_CER1-_AB", "cRH9_ACL_PD-_AB")))

# 使用 DotPlot 直接画图
dot_data <- DotPlot(ss22411, features = genes.to.plot, group.by = "group")$data

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

Idents(ss22411) <- ss22411$group
markers <- FindAllMarkers(ss22411, only.pos = T, logfc.threshold = 1)
write.table(markers,"ss22411markers.txt",sep = "\t")

VlnPlot(ss22411, group.by = "group",features = c("TTR","APOB"))
DotPlot(ss22411, features = "TTR", group.by = "group")$data

library(destiny)
DefaultAssay(ss22411)
mat <- as.matrix(ss22411@assays$RNA@data)
dms <- DiffusionMap(t(mat))

# 可视化：每个点颜色是 group / cell type
plot(dms$DC1, dms$DC2, col = ss22411$group, pch = 16,label = T)
legend("topright",                             # 图例位置
       legend = unique(ss22411$group),           # 图例标签
       col = unique(ss22411$group),              # 图例颜色（需与上面一致）
       pch = 16,                              # 图例图形
       title = "Group")



pseudobulk_expression_matrix2 <- Seurat::AverageExpression(ss22411,
                                                          group.by = "group", # Use the name of your metadata column
                                                          assays = "RNA",
                                                          slot = "counts" # CRITICAL: Get raw counts
)$RNA # Extract the matrix for the "RNA" assay


total_counts_per_sample2 <- colSums(pseudobulk_expression_matrix2)
pseudobulk_cpm2 <- sweep(pseudobulk_expression_matrix2,
                        MARGIN = 2,
                        STATS = total_counts_per_sample2 / 1000000,
                        FUN = "/")

filtered_pseudobulk_expression_matrix2 <- pseudobulk_cpm2[hvgs,]
correlation_matrix_cpm2 <- cor(filtered_pseudobulk_expression_matrix2, method = "spearman") # Or "spearman"
pheatmap(correlation_matrix_cpm2,
             display_numbers = TRUE,
             cluster_rows = FALSE, # 关闭行聚类
             cluster_cols = FALSE, # 关闭列聚类
             main = "Gene Expression Correlation (common genes, CPM)")

pdf("250610_correlationlss22411spearmanhgvs200022.pdf", width =7, height = 7)
pheatmap(correlation_matrix_cpm2,
         display_numbers = TRUE,
         cluster_rows = FALSE, # 关闭行聚类
         cluster_cols = FALSE, # 关闭列聚类
         color = my_colors_1,
         main = "Gene Expression Correlation ")
dev.off()

ss22411 <- FindVariableFeatures(ss22411, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(ss22411)

ss22411markers <- FindAllMarkers(ss22411, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
ss22411markers$ratio <- ss22411markers$pct.1/ss22411markers$pct.2
filtered_markersss22411 <- subset(ss22411markers, p_val_adj < 0.05)
new_markersss22411 <- add_cluster_cpm_to_markers(ss22411, filtered_markersss22411)
write.table(new_markersss22411,"250610markersss22411withCPM.txt",sep = "\t")
new_markersss22411<- read.table("250610markersss22411withCPM.txt",sep = "\t")

library(clusterProfiler)
library(org.Hs.eg.db)  # 人类基因注释，如果是小鼠用 org.Mm.eg.db
library(ggplot2)
library(dplyr)
library(stringr)
library(enrichplot)
library(patchwork)




plot_go_enrichment <- function(markers_df, 
                               species = "human",  # "human" 或 "mouse"
                               ont = "BP",         # "BP", "MF", "CC" 或 "ALL"
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               top_terms = 10,     # 每个cluster显示的top GO terms数量
                               plot_type = "combined", # "combined", "separate", "dotplot", "network"
                               min_genes = 5,      # 最少基因数过滤
                               logfc_cutoff = 0.25) {
  
  # 根据物种选择数据库
  if(species == "human") {
    orgdb <- org.Hs.eg.db
  } else if(species == "mouse") {
    orgdb <- org.Mm.eg.db
  } else {
    stop("species must be 'human' or 'mouse'")
  }
  
  # 过滤marker基因
  markers_filtered <- markers_df %>%
    filter(p_val_adj < pvalueCutoff, 
           abs(avg_log2FC) > logfc_cutoff,
           !is.na(gene))
  
  # 获取所有cluster
  clusters <- unique(markers_filtered$cluster)
  
  # 存储GO富集结果
  go_results <- list()
  
  # 对每个cluster进行GO富集分析
  for(cluster_id in clusters) {
    cat("Processing cluster:", cluster_id, "\n")
    
    # 获取该cluster的基因
    cluster_genes <- markers_filtered %>%
      filter(cluster == cluster_id) %>%
      pull(gene)
    
    if(length(cluster_genes) < min_genes) {
      cat("Cluster", cluster_id, "has fewer than", min_genes, "genes, skipping...\n")
      next
    }
    
    # 基因ID转换
    tryCatch({
      gene_ids <- bitr(cluster_genes, 
                       fromType = "SYMBOL",
                       toType = "ENTREZID", 
                       OrgDb = orgdb)
      
      if(nrow(gene_ids) == 0) {
        cat("No genes converted for cluster", cluster_id, "\n")
        next
      }
      
      # GO富集分析
      ego <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = orgdb,
                      ont = ont,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pvalueCutoff,
                      qvalueCutoff = qvalueCutoff,
                      readable = TRUE)
      
      if(!is.null(ego) && nrow(ego@result) > 0) {
        ego@result$cluster <- paste0("Cluster_", cluster_id)
        go_results[[paste0("cluster_", cluster_id)]] <- ego
      }
      
    }, error = function(e) {
      cat("Error processing cluster", cluster_id, ":", e$message, "\n")
    })
  }
  
  if(length(go_results) == 0) {
    stop("No GO enrichment results found for any cluster")
  }
  
  # 合并所有结果
  combined_results <- do.call(rbind, lapply(go_results, function(x) x@result))
  
  # 根据plot_type生成不同的图形
  if(plot_type == "combined") {
    # 组合柱状图
    plot_data <- combined_results %>%
      group_by(cluster) %>%
      slice_min(order_by = p.adjust, n = top_terms) %>%
      ungroup() %>%
      mutate(
        Description = str_wrap(Description, width = 50),
        log_pval = -log10(p.adjust),
        cluster = factor(cluster)
      )
    
    p <- ggplot(plot_data, aes(x = reorder(Description, log_pval), 
                               y = log_pval, 
                               fill = cluster)) +
      geom_col(position = "dodge", alpha = 0.8) +
      coord_flip() +
      facet_wrap(~cluster, scales = "free_y", ncol = 2) +
      labs(
        title = paste("GO", ont, "Enrichment Analysis by Cluster"),
        x = "GO Terms",
        y = "-log10(adjusted p-value)",
        fill = "Cluster"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none"
      ) +
      scale_fill_brewer(type = "qual", palette = "Set3")
    
  } else if(plot_type == "dotplot") {
    # 点图显示所有cluster
    plot_data <- combined_results %>%
      group_by(cluster) %>%
      slice_min(order_by = p.adjust, n = top_terms) %>%
      ungroup() %>%
      mutate(
        Description = str_wrap(Description, width = 40),
        GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
      )
    
    p <- ggplot(plot_data, aes(x = cluster, y = reorder(Description, -p.adjust))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
      scale_size_continuous(range = c(2, 8), name = "Gene Count") +
      labs(
        title = paste("GO", ont, "Enrichment Analysis - All Clusters"),
        x = "Cluster",
        y = "GO Terms"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8)
      )
    
  } else if(plot_type == "separate") {
    # 为每个cluster单独作图
    plot_list <- list()
    
    for(i in 1:length(go_results)) {
      cluster_name <- names(go_results)[i]
      ego <- go_results[[i]]
      
      if(nrow(ego@result) > 0) {
        p_cluster <- barplot(ego, showCategory = top_terms) +
          labs(title = paste("Cluster", gsub("cluster_", "", cluster_name))) +
          theme(axis.text.y = element_text(size = 8))
        
        plot_list[[cluster_name]] <- p_cluster
      }
      
    }
    
    # 使用patchwork合并图形
    if(length(plot_list) > 1) {
      p <- wrap_plots(plot_list, ncol = 2)
    } else {
      p <- plot_list[[1]]
    }
  }
  
  # 返回结果
  result <- list(
    plot = p,
    go_results = go_results,
    combined_results = combined_results,
    summary = combined_results %>%
      group_by(cluster) %>%
      summarise(
        n_terms = n(),
        top_term = Description[which.min(p.adjust)],
        min_pval = min(p.adjust),
        .groups = 'drop'
      )
  )
  
  return(result)
}

result <- plot_go_enrichment(new_markersss22411, 
                                                         species = "human",
                                                     ont = "BP",
                                                     plot_type = "dotplot")


emapplot(pairwise_termsim(result))                            
                             # # 显示图形
print(result$plot)
pp <- result$combined_results
write.table(pp, "250611GOresultdetails.txt",sep = "\t")
pdf("250610_ss22411GOterm.pdf", width =10, height = 7)
print(result$plot)
dev.off()

pca.stdev <- ss22411@reductions$pca@stdev  # 各PC的标准差
pca.var <- pca.stdev^2                 # 方差 = stdev^2
var.percent <- pca.var / sum(pca.var) * 100 
var.percent
