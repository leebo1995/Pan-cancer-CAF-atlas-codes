#### Fig 5A-C
g = qs::qread("CAF_and_Peri_harm.qs")
meta = qs::qread("edited_g_metadata.qs")
g = subset(g, cells = rownames(meta))
g@meta.data = meta

data = GetAssayData(g, assay = "RNA", slot = "counts")
cell_metadata = g@meta.data
gene_annotation = data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) = rownames(data)
cds = new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds = preprocess_cds(cds, num_dim = 30) 
cds.embed = cds@int_colData$reducedDims[[1]]
int.embed = Embeddings(g, reduction = "umap")
int.embed = int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims[["UMAP"]] = int.embed
cds = cluster_cells(cds, cluster_method = "louvain")
cds = learn_graph(cds) 
cds = order_cells(cds)
pdf("Monocle3.pdf", width = 10, height = 10)
plot_cells(cds, color_cells_by = "rename", cell_size = 0.5) + ggplot2::scale_color_manual(values = pal_rename)
plot_cells(cds, color_cells_by = "megacluster", cell_size = 0.5) + ggplot2::scale_color_manual(values = pal_megacluster)
plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1, label_groups_by_cluster = T)
dev.off()

AFD_genes <- c('PI16','COL1A1','HSPB6','HSPB8','CCL13','CXCL14','CSF2','CSF3','CXCL8','CCL19','CD74','CXCL10')
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$cell.type %in% c("AFD")]
AFD_lineage_cds <- order_cells(AFD_lineage_cds)

'CXCL2','HSPA1B','HSPA6','HSPH1','HSPA1A','HSP90AA1','HSPD1','HSP90AB1','HSPB6','HSPB8'
plot_genes_in_pseudotime(AFD_lineage_cds,color_cells_by="rename", min_expr=0.5, ncol= 2)+ ggplot2::scale_color_manual(values = colors_10）

plot_genes_in_pseudotime(cds[AFD_genes,])                                                                                                                     plot_genes_in_pseudotime(cds[AFD_genes,])

pdf("~/LB/data/pan cancer/CAF/fig/nsimCAF_monocle3/nsim_umap harmony.pdf", width = 11, height = 10)
DimPlot(g, label = T, group.by = "rename", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["rename"]][,1]))))) + ggtitle("0.1")
DimPlot(g, label = T, group.by = "rename1", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["rename1"]][,1]))))) + ggtitle("0.1")
DimPlot(g, label = T, group.by = "RNA_snn_res.0.1", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_snn_res.0.1"]][,1]))))) + ggtitle("0.1")
DimPlot(g, label = T, group.by = "RNA_snn_res.0.3", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_snn_res.0.3"]][,1]))))) + ggtitle("0.3")
DimPlot(g, label = T, group.by = "RNA_snn_res.0.4", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_snn_res.0.4"]][,1]))))) + ggtitle("0.4")
DimPlot(g, label = T, group.by = "RNA_snn_res.0.6", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_snn_res.0.6"]][,1]))))) + ggtitle("0.6")
DimPlot(g, label = T, group.by = "RNA_snn_res.0.8", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_snn_res.0.8"]][,1]))))) + ggtitle("0.8")
DimPlot(g, label = T, group.by = "RNA_snn_res.1.2", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_snn_res.1.2"]][,1]))))) + ggtitle("1.2")
DimPlot(g, label = F, group.by = "orig.ident", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["orig.ident"]][,1]))))) + NoLegend() + ggtitle("Batch")
DimPlot(g, label = F, group.by = "PatientID", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["PatientID"]][,1]))))) + NoLegend()
DimPlot(g, label = T, group.by = "Type", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["Type"]][,1])))))
DimPlot(g, label = T, group.by = "Tissue", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["Tissue"]][,1])))))
DimPlot(g, label = T, group.by = "GEO", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["GEO"]][,1])))))
DimPlot(g, label = T, group.by = "Gender", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["Gender"]][,1])))))
DimPlot(g, label = T, group.by = "Stage", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["Stage"]][,1])))))
DimPlot(g, label = T, group.by = "T", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["T"]][,1])))))
DimPlot(g, label = T, group.by = "N", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["N"]][,1])))))
DimPlot(g, label = T, group.by = "M", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["M"]][,1])))))
dev.off()



##Fig 5E 5G and Fig S15
library(limma)
library(edgeR)  # for voom function

# 假设counts已经加载
counts <- read.csv(file = '~/LB/data/RNA-seq/EN2_510_12h_24h_48h_72h_6d/data/counts_time.csv', header = TRUE, row.names = 1)

# 样本名称
sample <- c("NF1", "NF2", "NF3", "X12h1", "X12h2", "X12h3", "X24h1", "X24h2", "X24h3", "X48h1", "X48h2", "X48h3",
            "X72h1", "X72h2", "X72h3", "X6d1", "X6d2", "X6d3")

# 定义每个样本对应的条件（组）
group_labels <- c(
  rep("NF", 3), 
  rep("X12h", 3), 
  rep("X24h", 3), 
  rep("X48h", 3),
  rep("X72h", 3),
  rep("X6d", 3)
)

# 创建设计矩阵并转换为时间格式的列名
group <- factor(group_labels, levels = c("NF", "X12h", "X24h", "X48h", "X72h", "X6d"))
design <- model.matrix(~0 + group)

# 更改列名为时间点格式
colnames(design) <- c("time_0", "time_12", "time_24", "time_48", "time_72", "time_144")

# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 处理重复测量相关性
dupcor <- duplicateCorrelation(v, design, block=factor(sample))

# 如果duplicateCorrelation显示有显著的相关性，则在fit中包含这个参数
fit <- lmFit(v, design, block=factor(sample), correlation=dupcor$consensus)

# 对每个时间点进行对比分析
contrasts <- makeContrasts(
  'time_12_vs_time_0' = time_12 - time_0,
  'time_24_vs_time_0' = time_24 - time_0,
  'time_48_vs_time_0' = time_48 - time_0,
  'time_72_vs_time_0' = time_72 - time_0,
  'time_144_vs_time_0' = time_144 - time_0,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# 查找差异表达基因
de_genes <- topTable(fit2, adjust="fdr", number=Inf)

# 输出差异表达基因结果
write.csv(de_genes, file="differential_expression_results.csv")

# 构建对比矩阵以检测随时间的趋势
contrasts_trend <- makeContrasts(
  Time_Trend = time_12 + time_24 + time_48 + time_72 + time_144,
  levels = design
)

fit_trend <- contrasts.fit(fit, contrasts_trend)
fit_trend <- eBayes(fit_trend)

# 获取具有显著时间趋势的基因列表
trend_genes <- topTable(fit_trend, coef="Time_Trend", adjust="fdr", number=Inf)

# 输出时间趋势基因结果
write.csv(trend_genes, file="time_trend_genes.csv")

# 可视化结果
# 热图
heatmap.2(as.matrix(de_genes[order(abs(de_genes$logFC)), ]), col=heat.colors(256), scale="row", trace="none", main="Heatmap of Differentially Expressed Genes")

# 趋势图
plotMeans <- function(gene_id) {
  means <- aggregate(counts[gene_id,], by=list(group=sample_info$time), FUN=mean)
  matplot(means$group, t(means[, -1]), type='b', pch=19, col=rainbow(length(unique(means$group))), lty=1, xlab="Time (hours)", ylab="Expression Level", main=paste("Gene ID:", gene_id))
}

# 绘制前5个具有显著时间趋势的基因的趋势图
for (gene in rownames(trend_genes)[1:5]) {
  plotMeans(gene)
}


######多组差异通路图####
library(edgeR)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db) # 人类基因注释数据库
library(enrichplot)   # 用于绘制富集分析的结果

# 加载必要的包
library(GSVA)
library(limma)
library(msigdbr)

sample <- c("NF1", "NF2", "NF3", "X12h1", "X12h2", "X12h3", "X24h1", "X24h2", "X24h3", "X48h1", "X48h2", "X48h3",
            "X72h1", "X72h2", "X72h3", "X6d1", "X6d2", "X6d3")

# 定义每个样本对应的条件（组）
group_labels <- c(
  rep("NF", 3), 
  rep("X12h", 3), 
  rep("X24h", 3), 
  rep("X48h", 3),
  rep("X72h", 3),
  rep("X6d", 3)
)

# 创建设计矩阵并转换为时间格式的列名
group <- factor(group_labels, levels = c("NF", "X12h", "X24h", "X48h", "X72h", "X6d"))

# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge) # TMM 标准化
cpm_matrix <- cpm(dge)      # 计算每百万计数 (CPM)
expr_matrix <- log2(cpm_matrix + 1) # 对数转换

# 2. 下载通路基因集
msgdC2_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
msgdC5_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
msgdH <- msigdbr(species = "Homo sapiens", category = "H")

# 提取并过滤通路基因集
keggSet <- split(msgdC2_kegg$gene_symbol, msgdC2_kegg$gs_name)
gobpSet <- split(msgdC5_gobp$gene_symbol, msgdC5_gobp$gs_name)
hallmarkSet <- split(msgdH$gene_symbol, msgdH$gs_name)

common_genes <- intersect(rownames(expr_matrix), unlist(c(gobpSet, keggSet, hallmarkSet)))
expr_filtered <- expr_matrix[common_genes, ]

gobpSet_filtered <- lapply(gobpSet, function(gs) intersect(gs, common_genes))
keggSet_filtered <- lapply(keggSet, function(gs) intersect(gs, common_genes))
hallmarkSet_filtered <- lapply(hallmarkSet, function(gs) intersect(gs, common_genes))
pathwaySet <- c(keggSet_filtered, gobpSet_filtered, hallmarkSet_filtered)

# 3. GSVA 分析
zpar_pathwaySet <- zscoreParam(expr_filtered, geneSets = pathwaySet, minSize = 2)
gsva_results <- gsva(zpar_pathwaySet, verbose = TRUE)

# 定义一个函数来获取每个亚群的前12个通路
get_top12_per_cluster <- function(gsva_results, sample_info) {
  unique_clusters <- unique(sample_info$condition)
  top12_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, sample_info$condition == cluster]
    average_scores <- rowMeans(cluster_expr)
    top12 <- names(sort(average_scores, decreasing = TRUE))[1:12]
    top12_per_cluster[[cluster]] <- top12
  }
  
  return(top12_per_cluster)
}

# 获取每个亚群的前12个通路（所有通路一起考虑）
top12_per_cluster <- get_top12_per_cluster(gsva_results, sample_info)

# 构建新的数据矩阵，包含所有选定的通路
selected_pathways <- unique(unlist(top12_per_cluster))

# 选择对应的 GSVA 评分
selected_scores <- gsva_results[selected_pathways, ]

# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(selected_scores)) {
  stop("The number of samples does not match the number of columns in selected_scores.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("NF", "X12h", "X24h", "X48h", "X72h", "X6d")
sample_info_ordered <- sample_info[order(match(sample_info$condition, desired_order)), ]
selected_scores_ordered <- selected_scores[, colnames(selected_scores) %in% sample_info_ordered$sample]

# 聚合每个亚群的GSVA评分，并按照指定顺序排列
cluster_scores <- aggregate(t(selected_scores_ordered), by = list(cluster = sample_info_ordered$condition), FUN = mean)

# 按照指定顺序重新排序聚合结果
cluster_scores_ordered <- cluster_scores[match(desired_order, cluster_scores$cluster), ]

# 去掉聚合后的第一列（亚群名称）
heatmap_data <- as.matrix(cluster_scores_ordered[, -1])
rownames(heatmap_data) <- cluster_scores_ordered$cluster
colnames(heatmap_data) <- selected_pathways

# 绘制热图
pdf("~/LB/data/RNA-seq/EN2_510_12h_24h_48h_72h_6d/fig/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row", # 对每行进行标准化
  cluster_rows = FALSE,  # 关闭行聚类
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap (Top 12 Pathways per Cluster)",
  cellwidth = 20, # 设置每个单元格的宽度为60像素
  cellheight = 30 # 设置每个单元格的高度为40像素
)
dev.off()

# 将 gsva_results 转换为长格式，并添加组信息
gsva_long <- as.data.frame(t(gsva_results)) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "Pathway", values_to = "Score") %>%
  left_join(sample_info, by = c("Sample" = "sample"))

# 计算每个组的平均值
gsva_averages <- gsva_long %>%
  group_by(Pathway, condition) %>%
  summarise(AverageScore = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = AverageScore)


# 保存到 CSV 文件
output_path <- "~/LB/data/RNA-seq/EN2_510_12h_24h_48h_72h_6d/fig/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)


# 4. 定义要展示的39个通路
selected_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "GOBP_ORGANELLE_FISSION",
  "GOBP_CHROMOSOME_SEGREGATION",
  "GOBP_CELL_CYCLE_PROCESS",
  "GOBP_CELL_CYCLE",
  "HALLMARK_MYC_TARGETS_V1",
  "GOBP_RESPONSE_TO_HEAT",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
  "GOBP_CELLULAR_RESPONSE_TO_HEAT",
  "GOBP_RNA_PROCESSING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_MRNA_PROCESSING",
  "GOBP_MRNA_METABOLIC_PROCESS",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
  "GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",
  "GOBP_CELLULAR_RESPIRATION",
  "GOBP_CELLULAR_RESPONSE_TO_ENDOGENOUS_STIMULUS",
  "GOBP_RESPONSE_TO_ENDOGENOUS_STIMULUS",
  "GOBP_RESPONSE_TO_GROWTH_FACTOR",
  "GOBP_TISSUE_DEVELOPMENT",
  "GOBP_ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY",
  "GOBP_CELL_PART_MORPHOGENESIS",
  "GOBP_CELL_JUNCTION_ORGANIZATION",
  "GOBP_CELL_PROJECTION_ORGANIZATION",
  "GOBP_CELL_MORPHOGENESIS",
  "GOBP_REGULATION_OF_CELL_ADHESION",
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_CELL_MIGRATION",
  "GOBP_BIOLOGICAL_ADHESION",
  "GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES",
  "GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",
  "GOBP_RESPONSE_TO_OXYGEN_LEVELS",
  "GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "GOBP_CELL_CELL_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA"
)

# 筛选选定的通路
gsva_selected <- gsva_results[selected_pathways, ]

# 5. 构建新的数据矩阵，包含所有选定的通路
# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(gsva_selected)) {
  stop("The number of samples does not match the number of columns in gsva_selected.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("NF", "X12h", "X24h", "X48h", "X72h", "X6d")
sample_info_ordered <- sample_info[order(match(sample_info$condition, desired_order)), ]
gsva_selected_ordered <- gsva_selected[, colnames(gsva_selected) %in% sample_info_ordered$sample]

# 聚合每个亚群的GSVA评分，并按照指定顺序排列
cluster_scores <- aggregate(t(gsva_selected_ordered), by = list(cluster = sample_info_ordered$condition), FUN = mean)

# 按照指定顺序重新排序聚合结果
cluster_scores_ordered <- cluster_scores[match(desired_order, cluster_scores$cluster), ]

# 去掉聚合后的第一列（亚群名称）
heatmap_data <- as.matrix(cluster_scores_ordered[, -1])
rownames(heatmap_data) <- cluster_scores_ordered$cluster
colnames(heatmap_data) <- selected_pathways

# 绘制热图
pdf("~/LB/data/RNA-seq/EN2_510_12h_24h_48h_72h_6d/fig/GSVA/GSVA_heatmap_top39_pathways.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row", # 对每行进行标准化
  cluster_rows = FALSE,  # 关闭行聚类
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap (Top 39 Pathways)",
  cellwidth = 20, # 设置每个单元格的宽度为20像素
  cellheight = 30 # 设置每个单元格的高度为30像素
)
dev.off()










######BN_BCO_Time#####
######多组差异通路图####
library(edgeR)

library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db) # 人类基因注释数据库
library(enrichplot)   # 用于绘制富集分析的结果

# 加载必要的包
library(GSVA)
library(limma)
library(msigdbr)

counts <- read.csv(file = '~/LB/data/RNA-seq/BC/time/data/BN_counts_time.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型




sample <- c("BN1", "BN2", "BN3", "X12h1", "X12h2", "X12h3", "X24h1", "X24h2", "X24h3", "X48h1", "X48h2", "X48h3",
            "X72h1", "X72h2", "X72h3", "X6d1", "X6d2", "X6d3")

# 定义每个样本对应的条件（组）
group_labels <- c(
  rep("BN", 3), 
  rep("X12h", 3), 
  rep("X24h", 3), 
  rep("X48h", 3),
  rep("X72h", 3),
  rep("X6d", 3)
)

# 创建设计矩阵并转换为时间格式的列名
group <- factor(group_labels, levels = c("BN", "X12h", "X24h", "X48h", "X72h", "X6d"))

# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
counts[is.na(counts)] <- 0
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge) # TMM 标准化
cpm_matrix <- cpm(dge)      # 计算每百万计数 (CPM)
expr_matrix <- log2(cpm_matrix + 1) # 对数转换

# 2. 下载通路基因集
msgdC2_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
msgdC5_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
msgdH <- msigdbr(species = "Homo sapiens", category = "H")

# 提取并过滤通路基因集
keggSet <- split(msgdC2_kegg$gene_symbol, msgdC2_kegg$gs_name)
gobpSet <- split(msgdC5_gobp$gene_symbol, msgdC5_gobp$gs_name)
hallmarkSet <- split(msgdH$gene_symbol, msgdH$gs_name)

common_genes <- intersect(rownames(expr_matrix), unlist(c(gobpSet, keggSet, hallmarkSet)))
expr_filtered <- expr_matrix[common_genes, ]

gobpSet_filtered <- lapply(gobpSet, function(gs) intersect(gs, common_genes))
keggSet_filtered <- lapply(keggSet, function(gs) intersect(gs, common_genes))
hallmarkSet_filtered <- lapply(hallmarkSet, function(gs) intersect(gs, common_genes))
pathwaySet <- c(keggSet_filtered, gobpSet_filtered, hallmarkSet_filtered)

# 3. GSVA 分析
zpar_pathwaySet <- zscoreParam(expr_filtered, geneSets = pathwaySet, minSize = 2)
gsva_results <- gsva(zpar_pathwaySet, verbose = TRUE)

# 定义一个函数来获取每个亚群的前12个通路
get_top12_per_cluster <- function(gsva_results, sample_info) {
  unique_clusters <- unique(sample_info$condition)
  top12_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, sample_info$condition == cluster]
    average_scores <- rowMeans(cluster_expr)
    top12 <- names(sort(average_scores, decreasing = TRUE))[1:12]
    top12_per_cluster[[cluster]] <- top12
  }
  
  return(top12_per_cluster)
}

# 获取每个亚群的前12个通路（所有通路一起考虑）
top12_per_cluster <- get_top12_per_cluster(gsva_results, sample_info)

# 构建新的数据矩阵，包含所有选定的通路
selected_pathways <- unique(unlist(top12_per_cluster))

# 选择对应的 GSVA 评分
selected_scores <- gsva_results[selected_pathways, ]

# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(selected_scores)) {
  stop("The number of samples does not match the number of columns in selected_scores.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("BN", "X12h", "X24h", "X48h", "X72h", "X6d")
sample_info_ordered <- sample_info[order(match(sample_info$condition, desired_order)), ]
selected_scores_ordered <- selected_scores[, colnames(selected_scores) %in% sample_info_ordered$sample]

# 聚合每个亚群的GSVA评分，并按照指定顺序排列
cluster_scores <- aggregate(t(selected_scores_ordered), by = list(cluster = sample_info_ordered$condition), FUN = mean)

# 按照指定顺序重新排序聚合结果
cluster_scores_ordered <- cluster_scores[match(desired_order, cluster_scores$cluster), ]

# 去掉聚合后的第一列（亚群名称）
heatmap_data <- as.matrix(cluster_scores_ordered[, -1])
rownames(heatmap_data) <- cluster_scores_ordered$cluster
colnames(heatmap_data) <- selected_pathways

# 绘制热图
pdf("~/LB/data/RNA-seq/BC/time/fig/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row", # 对每行进行标准化
  cluster_rows = FALSE,  # 关闭行聚类
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap (Top 12 Pathways per Cluster)",
  cellwidth = 20, # 设置每个单元格的宽度为60像素
  cellheight = 30 # 设置每个单元格的高度为40像素
)
dev.off()

# 将 gsva_results 转换为长格式，并添加组信息
gsva_long <- as.data.frame(t(gsva_results)) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "Pathway", values_to = "Score") %>%
  left_join(sample_info, by = c("Sample" = "sample"))

# 计算每个组的平均值
gsva_averages <- gsva_long %>%
  group_by(Pathway, condition) %>%
  summarise(AverageScore = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = AverageScore)


# 保存到 CSV 文件
output_path <- "~/LB/data/RNA-seq/BC/time/fig/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)


# 4. 定义要展示的39个通路
selected_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "GOBP_ORGANELLE_FISSION",
  "GOBP_CHROMOSOME_SEGREGATION",
  "GOBP_CELL_CYCLE_PROCESS",
  "GOBP_CELL_CYCLE",
  "HALLMARK_MYC_TARGETS_V1",
  "GOBP_RESPONSE_TO_HEAT",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
  "GOBP_CELLULAR_RESPONSE_TO_HEAT",
  "GOBP_RNA_PROCESSING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_MRNA_PROCESSING",
  "GOBP_MRNA_METABOLIC_PROCESS",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
  "GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",
  "GOBP_CELLULAR_RESPIRATION",
  "GOBP_CELLULAR_RESPONSE_TO_ENDOGENOUS_STIMULUS",
  "GOBP_RESPONSE_TO_ENDOGENOUS_STIMULUS",
  "GOBP_RESPONSE_TO_GROWTH_FACTOR",
  "GOBP_TISSUE_DEVELOPMENT",
  "GOBP_ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY",
  "GOBP_CELL_PART_MORPHOGENESIS",
  "GOBP_CELL_JUNCTION_ORGANIZATION",
  "GOBP_CELL_PROJECTION_ORGANIZATION",
  "GOBP_CELL_MORPHOGENESIS",
  "GOBP_REGULATION_OF_CELL_ADHESION",
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_CELL_MIGRATION",
  "GOBP_BIOLOGICAL_ADHESION",
  "GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES",
  "GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",
  "GOBP_RESPONSE_TO_OXYGEN_LEVELS",
  "GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "GOBP_CELL_CELL_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA"
)

# 筛选选定的通路
gsva_selected <- gsva_results[selected_pathways, ]

# 5. 构建新的数据矩阵，包含所有选定的通路
# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(gsva_selected)) {
  stop("The number of samples does not match the number of columns in gsva_selected.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("BN", "X12h", "X24h", "X48h", "X72h", "X6d")
sample_info_ordered <- sample_info[order(match(sample_info$condition, desired_order)), ]
gsva_selected_ordered <- gsva_selected[, colnames(gsva_selected) %in% sample_info_ordered$sample]

# 聚合每个亚群的GSVA评分，并按照指定顺序排列
cluster_scores <- aggregate(t(gsva_selected_ordered), by = list(cluster = sample_info_ordered$condition), FUN = mean)

# 按照指定顺序重新排序聚合结果
cluster_scores_ordered <- cluster_scores[match(desired_order, cluster_scores$cluster), ]

# 去掉聚合后的第一列（亚群名称）
heatmap_data <- as.matrix(cluster_scores_ordered[, -1])
rownames(heatmap_data) <- cluster_scores_ordered$cluster
colnames(heatmap_data) <- selected_pathways

# 绘制热图
pdf("~/LB/data/RNA-seq/BC/time/fig/GSVA/GSVA_heatmap_top39_pathways.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row", # 对每行进行标准化
  cluster_rows = FALSE,  # 关闭行聚类
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap (Top 39 Pathways)",
  cellwidth = 20, # 设置每个单元格的宽度为20像素
  cellheight = 30 # 设置每个单元格的高度为30像素
)
dev.off()





######LN_LCO_Time#####
######多组差异通路图####
library(edgeR)

library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db) # 人类基因注释数据库
library(enrichplot)   # 用于绘制富集分析的结果

# 加载必要的包
library(GSVA)
library(limma)
library(msigdbr)

counts <- read.csv(file = '~/LB/data/RNA-seq/LC/time/data/LN_counts_time.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型




sample <- c("LN1", "LN2", "LN3", "X12h1", "X12h2", "X12h3", "X24h1", "X24h2", "X24h3", "X48h1", "X48h2", "X48h3",
            "X72h1", "X72h2", "X72h3", "X6d1", "X6d2", "X6d3")
group_labels <- c(rep("LN", 3), rep("X12h", 3), rep("X24h", 3), rep("X48h", 3), rep("X72h", 3), rep("X6d", 3))

# 创建设计矩阵并转换为时间格式的列名
group <- factor(group_labels, levels = c("LN", "X12h", "X24h", "X48h", "X72h", "X6d"))

# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
counts[is.na(counts)] <- 0
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge) # TMM 标准化
cpm_matrix <- cpm(dge)      # 计算每百万计数 (CPM)
expr_matrix <- log2(cpm_matrix + 1) # 对数转换

# 2. 下载通路基因集
msgdC2_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
msgdC5_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
msgdH <- msigdbr(species = "Homo sapiens", category = "H")

# 提取并过滤通路基因集
keggSet <- split(msgdC2_kegg$gene_symbol, msgdC2_kegg$gs_name)
gobpSet <- split(msgdC5_gobp$gene_symbol, msgdC5_gobp$gs_name)
hallmarkSet <- split(msgdH$gene_symbol, msgdH$gs_name)

common_genes <- intersect(rownames(expr_matrix), unlist(c(gobpSet, keggSet, hallmarkSet)))
expr_filtered <- expr_matrix[common_genes, ]

gobpSet_filtered <- lapply(gobpSet, function(gs) intersect(gs, common_genes))
keggSet_filtered <- lapply(keggSet, function(gs) intersect(gs, common_genes))
hallmarkSet_filtered <- lapply(hallmarkSet, function(gs) intersect(gs, common_genes))
pathwaySet <- c(keggSet_filtered, gobpSet_filtered, hallmarkSet_filtered)

# 3. GSVA 分析
zpar_pathwaySet <- zscoreParam(expr_filtered, geneSets = pathwaySet, minSize = 2)
gsva_results <- gsva(zpar_pathwaySet, verbose = TRUE)

# 定义一个函数来获取每个亚群的前12个通路
get_top12_per_cluster <- function(gsva_results, sample_info) {
  unique_clusters <- unique(sample_info$condition)
  top12_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, sample_info$condition == cluster]
    average_scores <- rowMeans(cluster_expr)
    top12 <- names(sort(average_scores, decreasing = TRUE))[1:12]
    top12_per_cluster[[cluster]] <- top12
  }
  
  return(top12_per_cluster)
}

# 获取每个亚群的前12个通路（所有通路一起考虑）
top12_per_cluster <- get_top12_per_cluster(gsva_results, sample_info)

# 构建新的数据矩阵，包含所有选定的通路
selected_pathways <- unique(unlist(top12_per_cluster))

# 选择对应的 GSVA 评分
selected_scores <- gsva_results[selected_pathways, ]

# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(selected_scores)) {
  stop("The number of samples does not match the number of columns in selected_scores.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("LN", "X12h", "X24h", "X48h", "X72h", "X6d")
sample_info_ordered <- sample_info[order(match(sample_info$condition, desired_order)), ]
selected_scores_ordered <- selected_scores[, colnames(selected_scores) %in% sample_info_ordered$sample]

# 聚合每个亚群的GSVA评分，并按照指定顺序排列
cluster_scores <- aggregate(t(selected_scores_ordered), by = list(cluster = sample_info_ordered$condition), FUN = mean)

# 按照指定顺序重新排序聚合结果
cluster_scores_ordered <- cluster_scores[match(desired_order, cluster_scores$cluster), ]

# 去掉聚合后的第一列（亚群名称）
heatmap_data <- as.matrix(cluster_scores_ordered[, -1])
rownames(heatmap_data) <- cluster_scores_ordered$cluster
colnames(heatmap_data) <- selected_pathways

# 绘制热图
pdf("~/LB/data/RNA-seq/LC/time/fig/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row", # 对每行进行标准化
  cluster_rows = FALSE,  # 关闭行聚类
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap (Top 12 Pathways per Cluster)",
  cellwidth = 20, # 设置每个单元格的宽度为60像素
  cellheight = 30 # 设置每个单元格的高度为40像素
)
dev.off()

# 将 gsva_results 转换为长格式，并添加组信息
gsva_long <- as.data.frame(t(gsva_results)) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "Pathway", values_to = "Score") %>%
  left_join(sample_info, by = c("Sample" = "sample"))

# 计算每个组的平均值
gsva_averages <- gsva_long %>%
  group_by(Pathway, condition) %>%
  summarise(AverageScore = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = AverageScore)


# 保存到 CSV 文件
output_path <- "~/LB/data/RNA-seq/LC/time/fig/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)


# 4. 定义要展示的39个通路
selected_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "GOBP_ORGANELLE_FISSION",
  "GOBP_CHROMOSOME_SEGREGATION",
  "GOBP_CELL_CYCLE_PROCESS",
  "GOBP_CELL_CYCLE",
  "HALLMARK_MYC_TARGETS_V1",
  "GOBP_RESPONSE_TO_HEAT",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
  "GOBP_CELLULAR_RESPONSE_TO_HEAT",
  "GOBP_RNA_PROCESSING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_MRNA_PROCESSING",
  "GOBP_MRNA_METABOLIC_PROCESS",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
  "GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",
  "GOBP_CELLULAR_RESPIRATION",
  "GOBP_CELLULAR_RESPONSE_TO_ENDOGENOUS_STIMULUS",
  "GOBP_RESPONSE_TO_ENDOGENOUS_STIMULUS",
  "GOBP_RESPONSE_TO_GROWTH_FACTOR",
  "GOBP_TISSUE_DEVELOPMENT",
  "GOBP_ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY",
  "GOBP_CELL_PART_MORPHOGENESIS",
  "GOBP_CELL_JUNCTION_ORGANIZATION",
  "GOBP_CELL_PROJECTION_ORGANIZATION",
  "GOBP_CELL_MORPHOGENESIS",
  "GOBP_REGULATION_OF_CELL_ADHESION",
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_CELL_MIGRATION",
  "GOBP_BIOLOGICAL_ADHESION",
  "GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES",
  "GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",
  "GOBP_RESPONSE_TO_OXYGEN_LEVELS",
  "GOBP_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "GOBP_CELL_CELL_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA"
)

# 筛选选定的通路
gsva_selected <- gsva_results[selected_pathways, ]

# 5. 构建新的数据矩阵，包含所有选定的通路
# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(gsva_selected)) {
  stop("The number of samples does not match the number of columns in gsva_selected.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("LN", "X12h", "X24h", "X48h", "X72h", "X6d")
sample_info_ordered <- sample_info[order(match(sample_info$condition, desired_order)), ]
gsva_selected_ordered <- gsva_selected[, colnames(gsva_selected) %in% sample_info_ordered$sample]

# 聚合每个亚群的GSVA评分，并按照指定顺序排列
cluster_scores <- aggregate(t(gsva_selected_ordered), by = list(cluster = sample_info_ordered$condition), FUN = mean)

# 按照指定顺序重新排序聚合结果
cluster_scores_ordered <- cluster_scores[match(desired_order, cluster_scores$cluster), ]

# 去掉聚合后的第一列（亚群名称）
heatmap_data <- as.matrix(cluster_scores_ordered[, -1])
rownames(heatmap_data) <- cluster_scores_ordered$cluster
colnames(heatmap_data) <- selected_pathways

# 绘制热图
pdf("~/LB/data/RNA-seq/LC/time/fig/GSVA/GSVA_heatmap_top39_pathways.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row", # 对每行进行标准化
  cluster_rows = FALSE,  # 关闭行聚类
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap (Top 39 Pathways)",
  cellwidth = 20, # 设置每个单元格的宽度为20像素
  cellheight = 30 # 设置每个单元格的高度为30像素
)
dev.off()







