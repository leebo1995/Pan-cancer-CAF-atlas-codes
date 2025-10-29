

######样本处理####
L1 <- Read10X(data.dir = "~/LB/data/mouse/data/L1/")
L2 <- Read10X(data.dir = "~/LB/data/mouse/data/L2/")
L3 <- Read10X(data.dir = "~/LB/data/mouse/data/L3/")
L4 <- Read10X(data.dir = "~/LB/data/mouse/data/L4/")

mouse_all <- readRDS(file = '~/LB/data/mouse/rds/mouse_all.rds')


samples <- list(
  Control = L1,Oxamate = L2,anti_PD1 = L3,Oxamate_anti_PD1 = L4)

# 创建 Seurat 对象列表
seurat_list <- lapply(names(samples), function(n) {
  obj <- CreateSeuratObject(counts = samples[[n]], project = n)
  obj$sample <- n
  obj$group <- n
  return(obj)
  })

combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(samples))

# 过滤
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)
# 归一化
combined <- NormalizeData(combined)

# 高变基因
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# 线性降维（可选整合）
combined <- ScaleData(combined)
combined <- RunPCA(combined)

library(harmony)
combined <- RunHarmony(combined, group.by.vars = 'orig.ident',assay.use = 'RNA', max.iter = 20)

combined <- FindNeighbors(combined, dims = 1:30, reduction = 'harmony')
combined <- FindClusters(combined, resolution = c(0.1,0.3,0.4,0.6,0.8,1.0,1.2))
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
DimPlot(combined, reduction='umap', label = TRUE, group.by = 'RNA_snn_res.1.2')

markers <- c("Cd3e", "Cd4", "Cd8a", "Cd79a", "Nkg7", "Ly6c2",
             "Adgre1", "Itgax", "S100a8", "Jchain", "Col1a1", "Pecam1", "Epcam","Krt7","Krt14")

# 绘制 FeaturePlot（多图排列）
FeaturePlot(combined, features = markers, cols = c("lightgrey", "red"), 
            reduction = "umap", ncol = 4)

combined@meta.data$group <-factor(combined@meta.data$group, levels = c('Control','Oxamate','anti_PD1','Oxamate_anti_PD1'))
DimPlot(combined, reduction = "umap", group.by = "group", split.by = "group")

DimPlot(combined, reduction='umap', label = TRUE, group.by = 'RNA_snn_res.0.4')

combined <- JoinLayers(combined)
Idents(combined) = "RNA_snn_res.0.4"
markers = FindAllMarkers(combined, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.7)
write.csv(markers,file = '~/LB/data/mouse/marker/markers_all_0.4.csv')

new_cluster_ids <- c(
  "0" = "Epithelial cells",      # Pcdh9, Slit2, Kirrel3
  "1" = "Epithelial cells",            # ENSMUSG 未知基因
  "2" = "Neutrophil",             # Cxcr2, G0s2
  "3" = "CD8+ T cells",              # Pdcd1, Klrh1
  "4" = "B cells",                # Pax5, Ms4a1
  "5" = "Macrophage",             # Mafb, Mertk
  "6" = "CD4+ T cells",           # Lef1, Tcf7
  "7" = "Epithelial cells", # Muc5ac, Pgc
  "8" = "Dendritic cells",        # Clec9a, Flt3
  "9" = "Macrophage",             # C1qa, C1qb
  "10" = "Acinar cells",          # Pnliprp1, Cpa1
  "11" = "NK cells",              # Ncr1, Klra8
  "12" = "Erythroid cells",       # Gypa, Hbb-bs
  "13" = "Neutrophil",            # Camp, Ngp
  "14" = "Fibroblast",            # Col5a2, Mfap5
  "15" = "Endothelial cells",     # Sox18, Gpihbp1
  "16" = "Plasma cells",          # Jchain, Igkc
  "17" = "Mast cells"             # Mrgprx2, Tpsab1
)

combined <- RenameIdents(combined, new_cluster_ids)
saveRDS (combined, file = '~/LB/data/mouse/rds/mouse_all.rds')

combined@meta.data$rename <- Idents(combined)
  
DimPlot(combined, reduction = "umap", group.by = "rename", split.by = "group")

#####计算Ro/e值###
install.packages("devtools")
devtools::install_github("Japrin/STARTRAC")
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)

data <- combined@meta.data
data <- data[, c("group", "rename")]
colnames(data) <- c("sample", "celltype")
data$tissue <- data$sample

Roe <- calTissueDist(data,
                     byPatient = FALSE,
                     colname.cluster = "celltype",
                     colname.patient = "sample",   # 这里是4个样本
                     colname.tissue = "tissue",    # 实际就是 sample
                     method = "chisq",
                     min.rowSum = 0)

Roe <- as.matrix(Roe)

breaks <- seq(0, max(Roe, na.rm = TRUE), length.out = 5)
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)),
                      c("#f6f8e6", "#f9a33e", "red"))

pdf("~/LB/data/mouse/fig/heatmap_all_Roe.pdf", width = 10, height = 10)
Heatmap(Roe,
        show_heatmap_legend = TRUE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        row_names_side = "right",
        column_names_side = "top",
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(
          title = "Ro/e value",
          at = breaks,
          labels = round(breaks, 2)
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y,
                    gp = gpar(fontsize = 10, col = "black"))
        })
dev.off()

pdf("~/LB/data/mouse/fig/mouse_all.pdf", width = 12, height = 10)
DimPlot(mouse_all, label = T, group.by = "rename", label.size = 4, pt.size = 1,  cols = pal_20)  + ggtitle("rename")
DimPlot(mouse_all, label = T, group.by = "group", label.size = 4, pt.size = 1,  cols = pal_20) + ggtitle("group")
dev.off()

pdf("~/LB/data/mouse/fig/mouse_all_split.pdf", width = 18, height = 10)
DimPlot(mouse_all, label = T, group.by = "rename", split.by = "group", label.size = 4, pt.size = 1,  cols = pal_20)+ ggtitle("split")
dev.off()

###提取巨噬细胞
mouse_Mac <- subset (mouse_all,  ident = c("Macrophage"))
mouse_Mac = NormalizeData(mouse_Mac)
mouse_Mac = FindVariableFeatures(mouse_Mac)
mouse_Mac = ScaleData(mouse_Mac)
mouse_Mac = RunPCA(mouse_Mac)
mouse_Mac <- RunHarmony(mouse_Mac, group.by.vars = 'orig.ident',assay.use = 'RNA', max.iter = 20)
mouse_Mac = FindNeighbors(mouse_Mac, dims = 1:30, reduction = "harmony")
mouse_Mac = FindClusters(mouse_Mac, resolution = c(0.1,0.3,0.4,0.6,0.8,1.2,1.5,1.8,2.0))
mouse_Mac = RunUMAP(mouse_Mac, reduction = "harmony", dims = 1:30, return.model = T)

pdf("~/LB/data/mouse/fig/mouse_all.pdf", width = 12, height = 10)
DimPlot(mouse_Mac, label = T, group.by = "RNA_snn_res.0.1", label.size = 4, pt.size = 1) 
DimPlot(mouse_Mac, label = T, group.by = "RNA_snn_res.0.3", label.size = 4, pt.size = 1)
DimPlot(mouse_Mac, label = T, group.by = "RNA_snn_res.0.4", label.size = 4, pt.size = 1) 
DimPlot(mouse_Mac, label = T, group.by = "RNA_snn_res.0.6", label.size = 4, pt.size = 1) 
DimPlot(mouse_Mac, label = T, group.by = "RNA_snn_res.0.8", label.size = 4, pt.size = 1) 
DimPlot(mouse_Mac, label = T, group.by = "RNA_snn_res.1.2", label.size = 4, pt.size = 1) 
DimPlot(mouse_Mac, label = T, group.by = "RNA_snn_res.1.5", label.size = 4, pt.size = 1) 
dev.off()




mouse_CAF1 <- JoinLayers(mouse_CAF1)
Idents(mouse_CAF1) = "RNA_snn_res.0.8"
markers = FindAllMarkers(mouse_CAF1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.7)



###提取CAF
mouse_CAF <- subset (combined,  ident = c("Fibroblast"))
mouse_CAF1 <- subset (combined,  ident = c("24"))

mouse_CAF1 = NormalizeData(mouse_CAF1)
mouse_CAF1 = FindVariableFeatures(mouse_CAF1)
mouse_CAF1 = ScaleData(mouse_CAF1)
mouse_CAF1 = RunPCA(mouse_CAF1)
mouse_CAF1 <- RunHarmony(mouse_CAF1, group.by.vars = 'orig.ident',assay.use = 'RNA', max.iter = 20)
mouse_CAF1 = FindNeighbors(mouse_CAF1, dims = 1:30, reduction = "harmony")
mouse_CAF1 = FindClusters(mouse_CAF1, resolution = c(0.1,0.3,0.4,0.6,0.8,1.2,1.5,1.8,2.0))
mouse_CAF1 = RunUMAP(mouse_CAF1, reduction = "harmony", dims = 1:30, return.model = T)
DimPlot(mouse_CAF1, reduction='umap', label = TRUE, group.by = 'RNA_snn_res.0.8')

mouse_CAF1 <- JoinLayers(mouse_CAF1)
Idents(mouse_CAF1) = "RNA_snn_res.0.8"
markers = FindAllMarkers(mouse_CAF1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.7)

write.csv(markers,file = '~/LB/data/mouse/marker/markers_mouse_CAF11_0.8.csv')


mouse_CAF1[["percent.mt"]] <- PercentageFeatureSet(mouse_CAF1, pattern = "^mt-")
VlnPlot(mouse_CAF1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mouse_CAF1 <- subset(mouse_CAF1, subset = nFeature_RNA < 5000 & percent.mt < 10)


mouse_CAF = NormalizeData(mouse_CAF)
mouse_CAF = FindVariableFeatures(mouse_CAF)
mouse_CAF = ScaleData(mouse_CAF)
mouse_CAF = RunPCA(mouse_CAF)
mouse_CAF <- RunHarmony(mouse_CAF, group.by.vars = 'orig.ident',assay.use = 'RNA', max.iter = 20)
mouse_CAF = FindNeighbors(mouse_CAF, dims = 1:30, reduction = "harmony")
mouse_CAF = FindClusters(mouse_CAF, resolution = c(0.1,0.3,0.4,0.6,0.8,1.2))
mouse_CAF = RunUMAP(mouse_CAF, reduction = "harmony", dims = 1:30, return.model = T)
DimPlot(mouse_CAF, reduction='umap', label = TRUE, group.by = 'RNA_snn_res.1.2')

mouse_CAF <- JoinLayers(mouse_CAF)
Idents(mouse_CAF) = "RNA_snn_res.1.2"
markers = FindAllMarkers(mouse_CAF, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.7)

write.csv(markers,file = '~/LB/data/mouse/marker/markers_mouse_CAF_1.21.csv')

mouse_CAF <- RenameIdents(mouse_CAF, `0`='lsCAF',`1`='Tem_iCAF/apCAF',`2`='Mac2_iCAF/Gran_iCAF',`3`='vCAF',`4`='iCAF',`5`='dCAF',`6`='myCAF',`7`='vCAF',`8`='dCAF')
mouse_CAF@meta.data$group<-factor(mouse_CAF@meta.data$group, levels = c("Control","Oxamate","anti_PD1","Oxamate_anti_PD1"))
mouse_CAF@meta.data$rename<-factor(mouse_CAF@meta.data$rename, levels = c("Mac2_iCAF/Gran_iCAF","Tem_iCAF/apCAF",'iCAF',"myCAF","dCAF","vCAF","lsCAF"))


mouse_CAF1 <- RenameIdents(mouse_CAF1, `0`='lsCAF',`1`='Tem_iCAF/apCAF',`2`='NF/Mono_iCAF',`3`='Mac2_iCAF/Gran_iCAF',`4`='Mac2_iCAF/Gran_iCAF',`5`='dCAF',`6`='myCAF',`7`='vCAF',`8`='dCAF')
mouse_CAF1@meta.data$rename <- Idents(mouse_CAF1)
mouse_CAF1@meta.data$rename<-factor(mouse_CAF1@meta.data$rename, levels = c("NF/Mono_iCAF","Mac2_iCAF/Gran_iCAF","Tem_iCAF/apCAF","myCAF","vCAF","dCAF","lsCAF"))
mouse_CAF1@meta.data$group<-factor(mouse_CAF1@meta.data$group, levels = c("Control","Oxamate","anti_PD1","Oxamate_anti_PD1"))

saveRDS (mouse_CAF1, file = '~/LB/data/mouse/rds/mouse_CAF1.rds')





DimPlot(mouse_CAF1, reduction = "umap", group.by = "rename", split.by = "group")


data <- mouse_CAF1@meta.data
data <- data[, c("group", "rename")]
colnames(data) <- c("sample", "celltype")
data$tissue <- data$sample

Roe <- calTissueDist(data,
                     byPatient = FALSE,
                     colname.cluster = "celltype",
                     colname.patient = "sample",   # 这里是4个样本
                     colname.tissue = "tissue",    # 实际就是 sample
                     method = "chisq",
                     min.rowSum = 0)

Roe <- as.matrix(Roe)

breaks <- seq(0, max(Roe, na.rm = TRUE), length.out = 5)
col_fun <- colorRamp2(c(min(Roe, na.rm = TRUE), 1, max(Roe, na.rm = TRUE)),
                      c("#f6f8e6", "#f9a33e", "red"))

pdf("~/LB/data/mouse/fig/heatmap_CAF_Roe1.pdf", width = 10, height = 10)
Heatmap(Roe,
        show_heatmap_legend = TRUE,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        row_names_side = "right",
        column_names_side = "top",
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(
          title = "Ro/e value",
          at = breaks,
          labels = round(breaks, 2)
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y,
                    gp = gpar(fontsize = 10, col = "black"))
        })
dev.off()

pdf("~/LB/data/mouse/fig/mouse_CAF1_1.pdf", width = 12, height = 10)
DimPlot(mouse_CAF1, label = T, group.by = "rename", label.size = 4, pt.size = 4,  cols = pal_20)  + ggtitle("rename")
DimPlot(mouse_CAF1, label = T, group.by = "group", label.size = 4, pt.size = 4,  cols = pal_20) + ggtitle("group")
dev.off()

pdf("~/LB/data/mouse/fig/mouse_CAF1_split_1.pdf", width = 18, height = 10)
DimPlot(mouse_CAF1, label = T, group.by = "rename", split.by = "group", label.size = 4, pt.size = 4,  cols = pal_20)+ ggtitle("split")
dev.off()


######多组差异通路图####
library(edgeR)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db) # 人类基因注释数据库
library(org.Mm.eg.db) # 人类基因注释数据库
library(enrichplot)   # 用于绘制富集分析的结果
library(GSVA)
library(limma)
library(msigdbr)

counts <- read.csv(file = '~/LB/data/mouse/RNA_seq/data/counts.csv', header = TRUE, row.names = 1)
sample <- c("Con1", "Con2", "Con3", "LDHA1", "LDHA2", "LDHA3", "PD1_1", "PD1_2", "PD1_3", "LDHA_PD1_1", "LDHA_PD1_2", "LDHA_PD1_3")

# 定义每个样本对应的条件（组）
group_labels <- c(
  rep("Con", 3), 
  rep("LDHA", 3), 
  rep("PD1_", 3), 
  rep("LDHA_PD1_", 3)
)

# 创建设计矩阵并转换为时间格式的列名
group <- factor(group_labels, levels = c("Con", "LDHA", "PD1_", "LDHA_PD1_"))

# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

#1. 读取并标准化表达数据 
# 假设 counts 是小鼠的 raw count matrix，行为 gene，列为样本
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
cpm_matrix <- cpm(dge)
expr_matrix <- log2(cpm_matrix + 1)

# ==== 2. 读取小鼠通路集 ====
msgdC2_kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
msgdC5_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
msgdH <- msigdbr(species = "Mus musculus", category = "H")

# 通路集转为 list 格式
keggSet <- split(msgdC2_kegg$gene_symbol, msgdC2_kegg$gs_name)
gobpSet <- split(msgdC5_gobp$gene_symbol, msgdC5_gobp$gs_name)
hallmarkSet <- split(msgdH$gene_symbol, msgdH$gs_name)

# 保留表达矩阵中存在的基因
common_genes <- intersect(rownames(expr_matrix), unlist(c(keggSet, gobpSet, hallmarkSet)))
expr_filtered <- expr_matrix[common_genes, ]

keggSet_filtered <- lapply(keggSet, function(x) intersect(x, common_genes))
gobpSet_filtered <- lapply(gobpSet, function(x) intersect(x, common_genes))
hallmarkSet_filtered <- lapply(hallmarkSet, function(x) intersect(x, common_genes))

pathwaySet <- c(keggSet_filtered, gobpSet_filtered, hallmarkSet_filtered)

# ==== 3. GSVA 分析 ====
zpar_pathwaySet <- zscoreParam(expr_filtered, geneSets = pathwaySet, minSize = 2)
gsva_results <- gsva(zpar_pathwaySet, verbose = TRUE)

# ==== 4. 获取每个组前12个通路 ====
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

top12_per_cluster <- get_top12_per_cluster(gsva_results, sample_info)
selected_pathways <- unique(unlist(top12_per_cluster))
selected_scores <- gsva_results[selected_pathways, ]

# ==== 5. 聚合按组绘制热图 ====
# 请提前设置 desired_order，比如：
desired_order <- c("Con", "LDHA", "PD1_", "LDHA_PD1_")

# 保证顺序一致
sample_info_ordered <- sample_info[order(match(sample_info$condition, desired_order)), ]
selected_scores_ordered <- selected_scores[, sample_info_ordered$sample]

cluster_scores <- aggregate(t(selected_scores_ordered), by = list(cluster = sample_info_ordered$condition), FUN = mean)
cluster_scores_ordered <- cluster_scores[match(desired_order, cluster_scores$cluster), ]
heatmap_data <- as.matrix(cluster_scores_ordered[, -1])
rownames(heatmap_data) <- cluster_scores_ordered$cluster
colnames(heatmap_data) <- selected_pathways

# ==== 6. 绘制热图 ====
pdf("~/LB/data/mouse/RNA_seq/fig/GSVA_heatmap_top12_per_cluster_mouse.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cellwidth = 20,
  cellheight = 30,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Heatmap - Top 12 Pathways per Mouse Cluster"
)
dev.off()

# ==== 7. 导出 CSV ====
gsva_long <- as.data.frame(t(gsva_results)) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "Pathway", values_to = "Score") %>%
  left_join(sample_info, by = c("Sample" = "sample"))

gsva_averages <- gsva_long %>%
  group_by(Pathway, condition) %>%
  summarise(AverageScore = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = AverageScore)

write.csv(gsva_averages, "~/LB/data/mouse/RNA_seq/fig/GSVA_scores_mouse_summary.csv", row.names = FALSE)


# 4. 定义要展示的39个通路
selected_pathways <- c(
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "HALLMARK_GLYCOLYSIS",
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
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
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
pdf("~/LB/data/mouse/RNA_seq/fig/GSVA_heatmap_top39_pathways.pdf", width = 20, height = 12)
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cellwidth = 20,
  cellheight = 30,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Heatmap - Top 12 Pathways per Mouse Cluster"
)
dev.off()



#Mac2_iCAF/Gran_iCAF各组的通路差异
library(Matrix)
library(GSVA)
library(limma)

# 1) 取 Mac2_iCAF/Gran_iCAF
mac2 <- subset(mouse_CAF1, idents = "Mac2_iCAF/Gran_iCAF")
exp <- GetAssayData(mac2, slot = "data")  # log-normalized
meta <- mac2@meta.data

# 2) 按 group 做伪Bulk（每组取行均值）
pb <- sapply(split(1:ncol(exp), meta$group), function(idx) Matrix::rowMeans(exp[, idx, drop=FALSE]))
# 3) GSVA/ssGSEA（推荐 ssGSEA 对样本数少时更稳）
gs_mouse <- msigdbr(species = "Mus musculus", category = "H")
gsets <- split(gs_mouse$gene_symbol, gs_mouse$gs_name)
paths <- c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
gsets <- gsets[paths]  # 同方案A
scores <- gsva(as.matrix(pb), gsets, method = "ssgsea", ssgsea.norm = TRUE)
scores <- zscoreParam(pb, geneSets = gsets, minSize = 2)
gsva_results <- gsva(scores, verbose = TRUE)

# 4) 可视化
pheatmap::pheatmap(scores, cluster_rows = TRUE, cluster_cols = FALSE,
                   main = "Mac2_iCAF/Gran_iCAF pathway (ssGSEA, pseudobulk)")
pheatmap(
  scores,
  scale = "row",
  cluster_rows = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cellwidth = 20,
  cellheight = 30,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Heatmap - Mac2_iCAF/Gran_iCAF per GROUP"
)

mat <- as.matrix(gsva_results)

## 2) （可选）按你四个GROUP顺序排列列
grp_order <- c("Control","Oxamate","anti_PD1","Oxamate_anti_PD1")
keep <- intersect(grp_order, colnames(mat))
mat  <- mat[, keep, drop = FALSE]

## 3) 行标准化（更好对比不同通路在各组的相对强弱）
mat_scaled <- t(scale(t(mat)))   # 每行减均值/除以sd；如不想标准化就用 mat

## 4) 画热图（行聚类、列不聚类；你也可以改成 cluster_cols=TRUE）

pdf("~/LB/data/mouse/fig/GSVA_Heatmap_Mac2iCAF_zscore_scaled.pdf", width=6, height=4)
pheatmap(mat_scaled,
         color = colorRampPalette(c("blue","white","red"))(50),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Mac2_iCAF/Gran_iCAF pathways (zscoreParam→GSVA)",
         fontsize_row = 10, fontsize_col = 11)
dev.off()
write.csv(mat_scaled, "~/LB/data/mouse/fig/GSVA_Mac2iCAF_zscore_scaled.csv")



## ========= 1) 取 Tem_iCAF/apCAF 子集 & 伪Bulk =========
tem <- subset(mouse_CAF1, idents = "Tem_iCAF/apCAF")
exp  <- GetAssayData(tem, slot = "data")              # genes x cells (log-normalized)
meta <- tem@meta.data

# 按组做伪Bulk（genes x groups）
pb <- sapply(split(1:ncol(exp), meta$group),
             function(idx) Matrix::rowMeans(exp[, idx, drop = FALSE]))
pb <- as.matrix(pb)

# 按你四组顺序排列列
grp_order <- c("Control","Oxamate","anti_PD1","Oxamate_anti_PD1")
keep <- intersect(grp_order, colnames(pb))
pb   <- pb[, keep, drop = FALSE]

## ========= 2) 准备免疫相关通路集合 =========
# Hallmark 免疫相关为主；另加“抗原呈递（MHC-II）”自定义集合
paths_immune <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"             # 可选：查看PD-1相关凋亡/活化信号
)

msig  <- msigdbr(species = "Mus musculus", category = "H")
gsets <- split(msig$gene_symbol, msig$gs_name)
gsets <- gsets[paths_immune]

# 自定义：抗原呈递/MHC-II（apCAF 特征）
custom_AP <- c("Cd74","H2-Aa","H2-Ab1","H2-Eb1","Ciita","H2-DMb1","Tap1","Tapbp","B2m","Ifitm3")
gsets$ANTIGEN_PRESENTATION_MHCII <- intersect(custom_AP, rownames(pb))

# 过滤各基因集，确保都在 pb 里
gsets <- lapply(gsets, function(v) intersect(v, rownames(pb)))
gsets <- gsets[ vapply(gsets, length, 1L) >= 2 ]       # 留至少2个基因的基因集

## ========= 3) 你的方法：zscoreParam → gsva() =========
# （提示：如果你只想跑 zscore 法，这里直接 method="zscore" 也可以）
scores_z <- zscoreParam(pb, geneSets = gsets, minSize = 2)  # 返回 GeneSetCollection z-scores
gsva_res <- gsva(scores_z, verbose = TRUE)                  # 通路 x 组 —— 这是要画图的矩阵
mat      <- as.matrix(gsva_res)

## ========= 4) 热图（行标准化 & 原始） =========
# 行标准化，便于比对各组相对强弱
mat_scaled <- t(scale(t(mat)))

pdf("~/LB/data/mouse/fig/GSVA_Heatmap_Tem_apCAF_scaled.pdf", width=6, height=4)
pheatmap(mat_scaled,
         color = colorRampPalette(c("blue","white","red"))(50),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Tem_apCAF pathways (zscoreParam→GSVA)",
         fontsize_row = 10, fontsize_col = 11)
dev.off()



write.csv(mat,        "~/LB/data/mouse/fig/GSVA_Tem_apCAF_raw.csv")
write.csv(mat_scaled, "~/LB/data/mouse/fig/GSVA_Tem_apCAF_scaled.csv")








