#######两两对比差异通路图####
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
library(ggplot2)


counts <- read.csv(file = '~/LB/data/RNA-seq/CSF2_CSF3/data/CSF2.csv', header = TRUE, row.names = 1)


# 样本名称和对应的条件
sample <- c("IgG2_1", "IgG2_2", "IgG2_3", "IgG2_4", "IgG2_5",
            "CSF2_1", "CSF2_2", "CSF2_3", "CSF2_4", "CSF2_5")
group_labels <- c(rep("IgG2", 5), rep("CSF2", 5))

# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
library(edgeR)
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

# 4. 创建设计矩阵，不含截距项
group <- factor(group_labels)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵
contrast.matrix <- makeContrasts(
  IgG2_vs_CSF2 = CSF2 - IgG2,
  levels = design
)

# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 7. 结果输出
results_list <- lapply(seq_along(colnames(contrast.matrix)), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})

names(results_list) <- colnames(contrast.matrix)

# 8. 可视化结果
df <- data.frame(ID = rownames(results_list[["IgG2_vs_CSF2"]]), score = results_list[["IgG2_vs_CSF2"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/CSF2_CSF3/fig/GSVA_IgG2vsCSF2.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "IgG2 vs CSF2") + # 确保标题是字符串
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

write.csv(results_list[["IgG2_vsCSF2_"]], file = "~/LB/data/RNA-seq/CSF2_CSF3/fig/GSVA_IgG2vsCSF2.csv", row.names = TRUE)



#######两两对比差异通路图####
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
library(ggplot2)


counts <- read.csv(file = '~/LB/data/RNA-seq/CSF2_CSF3/data/CSF3.csv', header = TRUE, row.names = 1)


# 样本名称和对应的条件
sample <- c("Con_1", "Con_2", "Con_3", "Con_4",
            "CSF3_1", "CSF3_2", "CSF3_3", "CSF3_4")
group_labels <- c(rep("Con", 4), rep("CSF3", 4))

# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
library(edgeR)
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

# 4. 创建设计矩阵，不含截距项
group <- factor(group_labels)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵
contrast.matrix <- makeContrasts(
  Con_vs_CSF3 = CSF3 - Con,
  levels = design
)

# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 7. 结果输出
results_list <- lapply(seq_along(colnames(contrast.matrix)), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})

names(results_list) <- colnames(contrast.matrix)

# 8. 可视化结果
df <- data.frame(ID = rownames(results_list[["Con_vs_CSF3"]]), score = results_list[["Con_vs_CSF3"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/CSF2_CSF3/fig/GSVA_ConvsCSF3.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "Con vs CSF3") + # 确保标题是字符串
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

write.csv(results_list[["Con_vs_CSF3"]], file = "~/LB/data/RNA-seq/CSF2_CSF3/fig/GSVA_ConvsCSF3.csv", row.names = TRUE)


