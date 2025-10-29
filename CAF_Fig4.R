
##Fig 4 and Fig S12 S13
######EN2_510_Mac#####
######两两对比火山图####
library(RColorBrewer)
rdylbu7 <-rev(brewer.pal(7,"RdYlBu"))



library(limma)
library(edgeR)  # for voom function

# 假设counts已经加载
counts <- read.csv(file = '~/LB/data/RNA-seq/EN2_510/data/counts_Mac.csv', header = TRUE, row.names = 1)

#样本名称
sample <- c("NF1", "NF2", "NF3", "NFMac1", "NFMac2", "NFMac3",
            "NF5101", "NF5102", "NF5103", "NF510Mac1", "NF510Mac2", "NF510Mac3")

# 定义每个样本对应的条件（组）
group_labels <- c(
  rep("NF", 3), # NF1, NF2, NF3
  rep("NFMac", 3), # NFMac1, NFMac2, NFMac3
  rep("NF510", 3), # NF5101, NF5102, NF5103
  rep("NF510Mac", 3) # NF510Mac1, NF510Mac2, NF510Mac3
)

# 创建设计矩阵
group <- factor(group_labels, levels = c("NF", "NFMac", "NF510", "NF510Mac"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 拟合线性模型
fit <- lmFit(v, design)

# 定义对比矩阵
contrasts_matrix <- makeContrasts(
  NF_vs_NF510 = NF510 - NF,
  NF_vs_NFMac = NFMac - NF,
  NF_vs_NF510Mac = NF510Mac - NF,
  NFMac_vs_NF510 = NF510 - NFMac,
  NF510_vs_NF510Mac = NF510Mac - NF510,
  NFMac_vs_NF510Mac = NF510Mac - NFMac,
  levels = design
)

# 应用对比
fit2 <- contrasts.fit(fit, contrasts_matrix)
eb <- eBayes(fit2)

# 检查对比矩阵是否成功应用并且有相应的系数
print(colnames(eb$coefficients))

# 创建一个函数来执行差异表达分析并绘制火山图
analyze_and_plot_limma <- function(contrast_name, eb) {
  # 确保对比名称存在于eb对象中
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  # 提取结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 如果结果为空，则跳过此对比
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  # 准备火山图的数据框
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  # 标记显著性以及增加颜色标记
  df$color <- "grey"  # 默认无显著性为灰色
  df$color[df$logFC > 1 ] <- "#E41A1C"  # 显著上调为红色
  df$color[df$logFC < -1 ] <- "#377EB8"  # 显著下调为蓝色
  
  # 提取logFC绝对值最大的前30位基因
  top_up <- head(df[order(-df$logFC), ], 30)  # logFC值最高的前30个基因
  top_down <- head(df[order(df$logFC), ], 30)  # logFC值最低的前30个基因
  
  # 创建一个包含所有要标注的基因的数据框
  top_genes <- rbind(top_up, top_down)
  
  # 火山图绘制
  p <-ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +  # 使用定义的颜色，不使用默认配色方案
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),  # 去掉主要网格线
      panel.grid.minor = element_blank(),  # 去掉次要网格线
      axis.line = element_line(color = "black"),  # 设置坐标轴线为黑色
      axis.ticks = element_line(color = "black")  # 设置刻度线颜色为黑色
      
    )
      +  # 移除图例
    
    # 使用geom_text_repel添加标签，确保标签不会重叠
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                             box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.5, "lines"),
                             segment.color = "grey50")
  
  return(p)
}

# 对每一对对比执行分析并绘制火山图
for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PDF文件路径
    pdf_path <- file.path("/home/data/t070206/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/Volcano/", paste0("Volcano Plot for ", contrast_name, ".pdf"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    
    # 开始绘制 PDF 文件
    pdf(pdf_path, width = 9, height = 6)
    
    # 打印图形到PDF
    print(plot)
    
    # 关闭PDF设备
    dev.off()
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}

for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PNG文件路径
    png_path <- file.path("/home/data/t070206/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/Volcano/", paste0("Volcano Plot for ", contrast_name, ".png"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
    
    # 使用 ggsave 保存 PNG 文件，指定宽度和高度以英寸为单位
    ggsave(filename = png_path, plot = plot, width = 9, height = 6, units = "in", dpi = 300)
    
    cat("对比", contrast_name, "已生成并保存为PNG图片\n")
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}



output_dir <- "/home/data/t070206/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/Volcano/"

# 确保输出目录存在，如果不存在则创建
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 对每一对比执行分析并导出为CSV文件
for (contrast_name in colnames(eb$coefficients)) {
  # 获取该对比的差异表达结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 构建完整的文件路径
  file_path <- file.path(output_dir, paste0(contrast_name, "_DEGs.csv"))
  
  # 将结果写入CSV文件
  write.csv(results, file=file_path, row.names=TRUE)
  
  cat("已保存文件:", file_path, "\n")
}




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

# 样本名称和对应的条件
sample <- c("NF1", "NF2", "NF3", "NFMac1", "NFMac2", "NFMac3",
            "NF5101", "NF5102", "NF5103", "NF510Mac1", "NF510Mac2", "NF510Mac3")
group_labels <- c(
  rep("NF", 3), # NF1, NF2,聂3
  rep("NFMac", 3), # NFMac1, NFMac2, NFMac3
  rep("NF510", 3), # NF5101, NF5102, NF5103
  rep("NF510Mac", 3) # NF510Mac1, NF510Mac2, NF510Mac3
)

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

#4.创建设计矩阵，不含截距项
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵 (替换 Mac 为 T)
contrast.matrix <- makeContrasts(
  NF_vs_NF510 = NF510 - NF,
  NF_vs_NFMac = NFMac - NF,
  NF_vs_NF510Mac = NF510Mac - NF,
  NFMac_vs_NF510 = NF510 - NFMac,
  NF510_vs_NF510Mac = NF510Mac - NF510,
  NFMac_vs_NF510Mac = NF510Mac - NFMac,
  levels = design
)

# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 7. 结果输出
results_list <- lapply(seq_along(contrast_names), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})

names(results_list) <- contrast_names

# 8. 可视化结果
  df <- data.frame(ID = rownames(results_list[["NF_vs_NFMac"]]), score = results_list[["NF_vs_NFMac"]]$logFC)
  df <- df[order(df$score), ]
  lowest_20 <- head(df, 20)
  highest_20 <- tail(df, 20)
  top_bottom_df <- rbind(lowest_20, highest_20)
  
  pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA_NFvsNFMac.pdf', width = 12, height = 8)
  ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
    labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFvsNFMac") + # 确保标题是字符串
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 0.6),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    coord_flip()
  dev.off()

  df <- data.frame(ID = rownames(results_list[["NF_vs_NF510"]]), score = results_list[["NF_vs_NF510"]]$logFC)
  df <- df[order(df$score), ]
  lowest_20 <- head(df, 20)
  highest_20 <- tail(df, 20)
  top_bottom_df <- rbind(lowest_20, highest_20)
  
  pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA_NFvsNF510.pdf', width = 12, height = 8)
  ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
    labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFvsNF510") + # 确保标题是字符串
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 0.6),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    coord_flip()
  dev.off()

  df <- data.frame(ID = rownames(results_list[["NF_vs_NF510Mac"]]), score = results_list[["NF_vs_NF510Mac"]]$logFC)
  df <- df[order(df$score), ]
  lowest_20 <- head(df, 20)
  highest_20 <- tail(df, 20)
  top_bottom_df <- rbind(lowest_20, highest_20)
  
  pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA_NFvsNF510Mac.pdf', width = 12, height = 8)
  ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
    labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFvsNF510Mac") + # 确保标题是字符串
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 0.6),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    coord_flip()
  dev.off()

  df <- data.frame(ID = rownames(results_list[["NFMac_vs_NF510"]]), score = results_list[["NFMac_vs_NF510"]]$logFC)
  df <- df[order(df$score), ]
  lowest_20 <- head(df, 20)
  highest_20 <- tail(df, 20)
  top_bottom_df <- rbind(lowest_20, highest_20)
  
  pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA_NFMacvsNF510.pdf', width = 12, height = 8)
  ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
    labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFMacvsNF510") + # 确保标题是字符串
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 0.6),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    coord_flip()
  dev.off()

  df <- data.frame(ID = rownames(results_list[["NFMac_vs_NF510Mac"]]), score = results_list[["NFMac_vs_NF510Mac"]]$logFC)
  df <- df[order(df$score), ]
  lowest_20 <- head(df, 20)
  highest_20 <- tail(df, 20)
  top_bottom_df <- rbind(lowest_20, highest_20)
  
  pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA_NFMacvsNF510Mac.pdf', width = 12, height = 8)
  ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
    labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFMacvsNF510Mac") + # 确保标题是字符串
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 0.6),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    coord_flip()
  dev.off()

  df <- data.frame(ID = rownames(results_list[["NF510_vs_NF510Mac"]]), score = results_list[["NF510_vs_NF510Mac"]]$logFC)
  df <- df[order(df$score), ]
  lowest_20 <- head(df, 20)
  highest_20 <- tail(df, 20)
  top_bottom_df <- rbind(lowest_20, highest_20)
  
  pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA_NF510vsNF510Mac.pdf', width = 12, height = 8)
  ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
    geom_bar(stat = "identity", alpha = 0.7) +
    scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
    labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NF510vsNF510Mac") + # 确保标题是字符串
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size = 0.6),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    coord_flip()
  dev.off()

write.csv(results_list[["NF_vs_NFMac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF_vs_NFMac.csv", row.names = TRUE)
write.csv(results_list[["NF_vs_NF510"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF_vs_NF510.csv", row.names = TRUE)
write.csv(results_list[["NF_vs_NF510Mac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF_vs_NF510Mac.csv", row.names = TRUE)
write.csv(results_list[["NFMac_vs_NF510"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NFMac_vs_NF510.csv", row.names = TRUE)
write.csv(results_list[["NFMac_vs_NF510Mac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NFMac_vs_NF510Mac.csv", row.names = TRUE)
write.csv(results_list[["NF510_vs_NF510Mac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF510_vs_NF510Mac.csv", row.names = TRUE)



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

# 样本名称和对应的条件 (替换 Mac 为 T)
sample <- c("NF1", "NF2", "NF3", "NFMac1", "NFMac2", "NFMac3",
            "NF5101", "NF5102", "NF5103", "NF510Mac1", "NF510Mac2", "NF510Mac3")
group_labels <- c(
  rep("NF", 3), # NF1, NF2, NF3
  rep("NFMac", 3), # NFMac1, NFMac2, NFMac3
  rep("NF510", 3), # NF5101, NF5102, NF5103
  rep("NF510Mac", 3) # NF510Mac1, NF510Mac2, NF510Mac3
)

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
desired_order <- c("NF", "NFMac", "NF510", "NF510Mac")
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
pdf("~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 14, height = 12)
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
output_path <- "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)
######EN2_510_T#####
#######两两对比火山图#####
library(RColorBrewer)
rdylbu7 <- rev(brewer.pal(7, "RdYlBu"))

library(limma)
library(edgeR)  # for voom function

# 假设counts已经加载
counts <- read.csv(file = '~/LB/data/RNA-seq/EN2_510/data/counts_T.csv', header = TRUE, row.names = 1)

# 样本名称 (替换 Mac 为 T)
sample <- c("NF1", "NF2", "NF3", "NFT1", "NFT2", "NFT3",
            "NF5101", "NF5102", "NF5103", "NF510T1", "NF510T2", "NF510T3")

# 定义每个样本对应的条件（组）(替换 Mac 为 T)
group_labels <- c(
  rep("NF", 3), # NF1, NF2, NF3
  rep("NFT", 3), # NFT1, NFT2, NFT3
  rep("NF510", 3), # NF5101, NF5102, NF5103
  rep("NF510T", 3) # NF510T1, NF510T2, NF510T3
)

# 创建设计矩阵
group <- factor(group_labels, levels = c("NF", "NFT", "NF510", "NF510T"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 拟合线性模型
fit <- lmFit(v, design)

# 定义对比矩阵 (替换 Mac 为 T)
contrasts_matrix <- makeContrasts(
  NF_vs_NF510 = NF510 - NF,
  NF_vs_NFT = NFT - NF,
  NF_vs_NF510T = NF510T - NF,
  NFT_vs_NF510 = NF510 - NFT,
  NF510_vs_NF510T = NF510T - NF510,
  NFT_vs_NF510T = NF510T - NFT,
  levels = design
)

# 应用对比
fit2 <- contrasts.fit(fit, contrasts_matrix)
eb <- eBayes(fit2)



# 函数 analyze_and_plot_limma 保持不变，因为它不直接涉及 "Mac"
analyze_and_plot_limma <- function(contrast_name, eb) {
  # 确保对比名称存在于eb对象中
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  # 提取结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 如果结果为空，则跳过此对比
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  # 准备火山图的数据框
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  # 标记显著性以及增加颜色标记
  df$color <- "grey"  # 默认无显著性为灰色
  df$color[df$logFC > 1 ] <- "#E41A1C"  # 显著上调为红色
  df$color[df$logFC < -1 ] <- "#377EB8"  # 显著下调为蓝色
  
  # 提取logFC绝对值最大的前30位基因
  top_up <- head(df[order(-df$logFC), ], 30)  # logFC值最高的前30个基因
  top_down <- head(df[order(df$logFC), ], 30)  # logFC值最低的前30个基因
  
  # 创建一个包含所有要标注的基因的数据框
  top_genes <- rbind(top_up, top_down)
  
  # 火山图绘制
  p <-ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +  # 使用定义的颜色，不使用默认配色方案
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),  # 去掉主要网格线
      panel.grid.minor = element_blank(),  # 去掉次要网格线
      axis.line = element_line(color = "black"),  # 设置坐标轴线为黑色
      axis.ticks = element_line(color = "black")  # 设置刻度线颜色为黑色
      
    )
  +  # 移除图例
    
    # 使用geom_text_repel添加标签，确保标签不会重叠
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                             box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.5, "lines"),
                             segment.color = "grey50")
  
  return(p)
}




# 对每一对对比执行分析并绘制火山图 (使用新的对比名称)
for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PDF文件路径
    pdf_path <- file.path("/home/data/t070206/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/Volcano/", paste0("Volcano Plot for ", contrast_name, ".pdf"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    
    # 开始绘制 PDF 文件
    pdf(pdf_path, width = 9, height = 6)
    
    # 打印图形到PDF
    print(plot)
    
    # 关闭PDF设备
    dev.off()
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}

for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PNG文件路径
    png_path <- file.path("/home/data/t070206/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/Volcano/", paste0("Volcano Plot for ", contrast_name, ".png"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
    
    # 使用 ggsave 保存 PNG 文件，指定宽度和高度以英寸为单位
    ggsave(filename = png_path, plot = plot, width = 9, height = 6, units = "in", dpi = 300)
    
    cat("对比", contrast_name, "已生成并保存为PNG图片\n")
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}




output_dir <- "/home/data/t070206/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/Volcano/"

# 确保输出目录存在，如果不存在则创建
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 对每一对比执行分析并导出为CSV文件
for (contrast_name in colnames(eb$coefficients)) {
  # 获取该对比的差异表达结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 构建完整的文件路径
  file_path <- file.path(output_dir, paste0(contrast_name, "_DEGs.csv"))
  
  # 将结果写入CSV文件
  write.csv(results, file=file_path, row.names=TRUE)
  
  cat("已保存文件:", file_path, "\n")
}





#######两两对比差异通路图#######
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

# 样本名称和对应的条件 (替换 Mac 为 T)
sample <- c("NF1", "NF2", "NF3", "NFT1", "NFT2", "NFT3",
            "NF5101", "NF5102", "NF5103", "NF510T1", "NF510T2", "NF510T3")
group_labels <- c(
  rep("NF", 3), # NF1, NF2, NF3
  rep("NFT", 3), # NFT1, NFT2, NFT3
  rep("NF510", 3), # NF5101, NF5102, NF5103
  rep("NF510T", 3) # NF510T1, NF510T2, NF510T3
)

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

# 4. 构建设计矩阵
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵 (替换 Mac 为 T)
contrast.matrix <- makeContrasts(
  NF_vs_NF510 = NF510 - NF,
  NF_vs_NFT = NFT - NF,
  NF_vs_NF510T = NF510T - NF,
  NFT_vs_NF510 = NF510 - NFT,
  NF510_vs_NF510T = NF510T - NF510,
  NFT_vs_NF510T = NF510T - NFT,
  levels = design
)

# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 获取对比名称列表
contrast_names <- colnames(coef(fit2))

# 7. 结果输出
results_list <- lapply(seq_along(contrast_names), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})

names(results_list) <- contrast_names

# 8. 可视化结果

df <- data.frame(ID = rownames(results_list[["NF_vs_NFT"]]), score = results_list[["NF_vs_NFT"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_NFvsNFT.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFvsNFT") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

df <- data.frame(ID = rownames(results_list[["NF_vs_NF510"]]), score = results_list[["NF_vs_NF510"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_NFvsNF510.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFvsNF510") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

df <- data.frame(ID = rownames(results_list[["NF_vs_NF510T"]]), score = results_list[["NF_vs_NF510T"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_NFvsNF510T.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFvsNF510T") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

df <- data.frame(ID = rownames(results_list[["NFT_vs_NF510"]]), score = results_list[["NFT_vs_NF510"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_NFTvsNF510.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFTvsNF510") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

df <- data.frame(ID = rownames(results_list[["NF510_vs_NF510T"]]), score = results_list[["NF510_vs_NF510T"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_NF510vsNF510T.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NF510vsNF510T") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

df <- data.frame(ID = rownames(results_list[["NFT_vs_NF510T"]]), score = results_list[["NFT_vs_NF510T"]]$logFC)
df <- df[order(df$score), ]
lowest_20 <- head(df, 20)
highest_20 <- tail(df, 20)
top_bottom_df <- rbind(lowest_20, highest_20)

pdf('~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_NFTvsNF510T.pdf', width = 12, height = 8)
ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
  geom_bar(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
  labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = "NFTvsNF510T") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_flip()
dev.off()

# 保存差异分析结果到CSV文件
write.csv(results_list[["NF_vs_NF510"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/NF_vs_NF510.csv", row.names = TRUE)
write.csv(results_list[["NF_vs_NFT"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/NF_vs_NFT.csv", row.names = TRUE)
write.csv(results_list[["NF_vs_NF510T"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/NF_vs_NF510T.csv", row.names = TRUE)
write.csv(results_list[["NFT_vs_NF510"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/NFT_vs_NF510.csv", row.names = TRUE)
write.csv(results_list[["NF510_vs_NF510T"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/NF510_vs_NF510T.csv", row.names = TRUE)
write.csv(results_list[["NFT_vs_NF510T"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/NFT_vs_NF510T.csv", row.names = TRUE)




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

# 样本名称和对应的条件 (替换 Mac 为 T)
sample <- c("NF1", "NF2", "NF3", "NFT1", "NFT2", "NFT3",
            "NF5101", "NF5102", "NF5103", "NF510T1", "NF510T2", "NF510T3")
group_labels <- c(
  rep("NF", 3), # NF1, NF2, NF3
  rep("NFT", 3), # NFT1, NFT2, NFT3
  rep("NF510", 3), # NF5101, NF5102, NF5103
  rep("NF510T", 3) # NF510T1, NF510T2, NF510T3
)

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
desired_order <- c("NF", "NFT", "NF510", "NF510T")
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
pdf("~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 14, height = 12)
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
output_path <- "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_T/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)


#转录组热图
markers <- c('PI16','CD34','ADH1B','THBS4','CLU','CCDC80','COL10A1','COL11A1','COL1A1','MMP11','MMP1','MMP3','CSF3','CXCL3','CXCL5','IL1B','IL1A','IL6','CXCL8','CSF2','CXCL14','CCL2','CCL8','CCL11','CCL13','CXCL9','CXCL10','CXCL11','CCL19','CCL21','CD74','HLA−DQA1','HLA−DPA1')

BN_BC_Mac <- read.csv(file = '~/LB/data/RNA-seq/BC/co/data/FPKM_Mac.csv', header = TRUE, row.names = 1)
BN_BC_Mac=as.matrix(BN_BC_Mac)  


markers <- gsub("−", "-", markers)
row_Mac <- BN_BC_Mac[markers, ]
row_Mac=as.matrix(row_Mac)  

pheatmap(row_mac, scale = "row",cluster_rows = F,cluster_cols = F,cellwidth = 40,cellheight = 20,border='black',show_rownames=TRUE,show_colnames=TRUE,color = rdylbu7)


# 读取数据
BN_BC_Mac <- read.csv(file = "~/LB/data/RNA-seq/BC/co/data/FPKM_Mac_1.csv", 
                      header = TRUE, 
                      row.names = 1)
BN_BC_Mac <- as.matrix(BN_BC_Mac)

# 定义markers
markers <- c('PI16','CD34','ADH1B','THBS4','CLU','CCDC80','COL10A1','COL11A1',
             'COL1A1','MMP11','MMP1','MMP3','CSF3','CXCL3','CXCL5','IL1B',
             'IL1A','IL6','CXCL8','CSF2','CXCL14','CCL2','CCL8','CCL11',
             'CCL13','CXCL9','CXCL10','CXCL11','CCL19','CCL21','CD74')

# 提取指定行
row_Mac <- BN_BC_Mac[markers, ]

# 将字符矩阵转换为数值矩阵
row_Mac_numeric <- type.convert(row_Mac, as.is = TRUE)
row_Mac_numeric
# 绘制热图
pdf("~/LB/data/RNA-seq/BC/co/fig/BN_BC_mac_heatmap.pdf", width = 10, height = 12)
pheatmap(row_Mac_numeric, 
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         cellheight = 20,
         border = "black",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rdylbu7)
dev.off()

pdf("~/LB/data/RNA-seq/BC/co/fig/BN_BC_mac_heatmap.pdf", width = 10, height = 12)
pheatmap(BN_BC_Mac, 
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         cellheight = 20,
         border = "black",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rdylbu7)
dev.off()


BN_BC_T <- read.csv(file = "~/LB/data/RNA-seq/BC/co/data/FPKM_T_1.csv", 
                      header = TRUE, 
                      row.names = 1)
BN_BC_T <- as.matrix(BN_BC_T)

pdf("~/LB/data/RNA-seq/BC/co/fig/BN_BC_T_heatmap.pdf", width = 10, height = 12)
pheatmap(BN_BC_T, 
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         cellheight = 20,
         border = "black",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rdylbu7)
dev.off()


BN_BC_time <- read.csv(file = "~/LB/data/RNA-seq/BC/time/data/BC_FPKM-time_1.csv", 
                    header = TRUE, 
                    row.names = 1)
BN_BC_time <- as.matrix(BN_BC_time)

pdf("~/LB/data/RNA-seq/BC/time/fig/BN_BC_time_heatmap1.pdf", width = 10, height = 12)
pheatmap(BN_BC_time, 
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         cellheight = 20,
         border = "black",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rdylbu11)
dev.off()



LN_LC_Mac <- read.csv(file = "~/LB/data/RNA-seq/LC/co/data/FKPM_Mac.csv", 
                      header = TRUE, 
                      row.names = 1)
LN_LC_Mac <- as.matrix(LN_LC_Mac)

pdf("~/LB/data/RNA-seq/LC/co/fig/LN_LC_mac_heatmap.pdf", width = 10, height = 12)
pheatmap(LN_LC_Mac, 
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         cellheight = 20,
         border = "black",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rdylbu7)
dev.off()


LN_LC_T <- read.csv(file = "~/LB/data/RNA-seq/LC/co/data/FKPM_T.csv", 
                      header = TRUE, 
                      row.names = 1)
LN_LC_T <- as.matrix(LN_LC_T)

pdf("~/LB/data/RNA-seq/LC/co/fig/LN_LC_T_heatmap.pdf", width = 10, height = 12)
pheatmap(LN_LC_T, 
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         cellheight = 20,
         border = "black",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rdylbu7)
dev.off()

LN_LC_time <- read.csv(file = "~/LB/data/RNA-seq/LC/time/data/LC_FPKM_time_2.csv", 
                       header = TRUE, 
                       row.names = 1)
LN_LC_time <- as.matrix(LN_LC_time)

pdf("~/LB/data/RNA-seq/LC/time/fig/LN_LC_time_heatmap1.pdf", width = 10, height = 12)
pheatmap(LN_LC_time, 
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         cellheight = 20,
         border = "black",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = rdylbu11)
dev.off()










######BN_BCO_Mac#####
######两两对比火山图####
library(RColorBrewer)
rdylbu7 <-rev(brewer.pal(7,"RdYlBu"))



library(limma)
library(edgeR)  # for voom function

# 假设counts已经加载
counts <- read.csv(file = '~/LB/data/RNA-seq/BC/co/data/BN_counts_Mac.csv', header = TRUE, row.names = 1)

# 样本名称修改（NF→BN，510→BCO）
sample <- c("BN1", "BN2", "BN3", "BNMac1", "BNMac2", "BNMac3",
            "BNBCO1", "BNBCO2", "BNBCO3", "BNBCOMac1", "BNBCOMac2", "BNBCOMac3")

# 组标签修改 
group_labels <- c(
  rep("BN", 3),        # 原NF组 
  rep("BNMac", 3),     # 原NFMac组 
  rep("BNBCO", 3),     # 原NF510组 
  rep("BNBCOMac", 3)   # 原NF510Mac组 
)

# 因子水平修改 
group <- factor(group_labels, levels = c("BN", "BNMac", "BNBCO", "BNBCOMac"))


design <- model.matrix(~0 + group)
colnames(design) <- levels(group)



counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型
counts[is.na(counts)] <- 0
v <- voom(counts, design, plot = FALSE)

# 拟合线性模型
fit <- lmFit(v, design)

# 对比矩阵修改 
contrasts_matrix <- makeContrasts(
  BN_vs_BNBCO = BNBCO - BN,
  BN_vs_BNMac = BNMac - BN,
  BN_vs_BNBCOMac = BNBCOMac - BN,
  BNMac_vs_BNBCO = BNBCO - BNMac,
  BNBCO_vs_BNBCOMac = BNBCOMac - BNBCO,
  BNMac_vs_BNBCOMac = BNBCOMac - BNMac,
  levels = design 
)
# 应用对比
fit2 <- contrasts.fit(fit, contrasts_matrix)
eb <- eBayes(fit2)

# 检查对比矩阵是否成功应用并且有相应的系数
print(colnames(eb$coefficients))

# 创建一个函数来执行差异表达分析并绘制火山图

analyze_and_plot_limma <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                             box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.5, "lines"),
                             segment.color = "grey50")
  
  return(p)
}

# 对每一对对比执行分析并绘制火山图
for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PDF文件路径
    pdf_path <- file.path("/home/data/t070206/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/Volcano/", paste0("Volcano Plot for ", contrast_name, ".pdf"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    
    # 开始绘制 PDF 文件
    pdf(pdf_path, width = 9, height = 6)
    
    # 打印图形到PDF
    print(plot)
    
    # 关闭PDF设备
    dev.off()
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}


analyze_and_plot_limma1 <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) 
  
  return(p)
}


for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma1(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PNG文件路径
    png_path <- file.path("/home/data/t070206/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/Volcano/", paste0("Volcano Plot for ", contrast_name, ".png"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
    
    # 使用 ggsave 保存 PNG 文件，指定宽度和高度以英寸为单位
    ggsave(filename = png_path, plot = plot, width = 9, height = 6, units = "in", dpi = 300)
    
    cat("对比", contrast_name, "已生成并保存为PNG图片\n")
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}



output_dir <- "/home/data/t070206/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/Volcano/"

# 确保输出目录存在，如果不存在则创建
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 对每一对比执行分析并导出为CSV文件
for (contrast_name in colnames(eb$coefficients)) {
  # 获取该对比的差异表达结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 构建完整的文件路径
  file_path <- file.path(output_dir, paste0(contrast_name, "_DEGs.csv"))
  
  # 将结果写入CSV文件
  write.csv(results, file=file_path, row.names=TRUE)
  
  cat("已保存文件:", file_path, "\n")
}






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

counts <- read.csv(file = '~/LB/data/RNA-seq/BC/co/data/BN_counts_Mac.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型


# 样本名称和对应的条件
sample <- c("BN1", "BN2", "BN3", "BNMac1", "BNMac2", "BNMac3",
            "BNBCO1", "BNBCO2", "BNBCO3", "BNBCOMac1", "BNBCOMac2", "BNBCOMac3")

# 组标签修改 
group_labels <- c(
  rep("BN", 3),        # 原NF组 
  rep("BNMac", 3),     # 原NFMac组 
  rep("BNBCO", 3),     # 原NF510组 
  rep("BNBCOMac", 3)   # 原NF510Mac组 
)

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

#4.创建设计矩阵，不含截距项
group <- factor(group_labels, levels = c("BN", "BNMac", "BNBCO", "BNBCOMac"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵 (替换 Mac 为 T)
contrast.matrix <- makeContrasts(
  BNvsBNBCO = BNBCO - BN,
  BNvsBNMac = BNMac - BN,
  BNvsBNBCOMac = BNBCOMac - BN,
  BNMacvsBNBCO = BNBCO - BNMac,
  BNBCOvsBNBCOMac = BNBCOMac - BNBCO,
  BNMacvsBNBCOMac = BNBCOMac - BNMac,
  levels = design
)

# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#7. 结果输出
contrast_names <- colnames(contrast.matrix)
contrast_list <- contrast_names  # 用于后续循环
results_list <- lapply(seq_along(contrast_names), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})
names(results_list) <- contrast_names

fig_dir <- "~/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/GSVA/"
csv_dir <- file.path(fig_dir, "GSVA")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 循环处理每组
for (contrast in contrast_list) {
  df <- data.frame(ID = rownames(results_list[[contrast]]), score = results_list[[contrast]]$logFC)
  df <- df[order(df$score), ]
  lowest_10 <- head(df, 10)
  highest_10 <- tail(df, 10)
  top_bottom_df <- rbind(lowest_10, highest_10)  # ✅ 修复这里
  
  # 画图
  pdf(file.path(fig_dir, paste0("GSVA_", contrast, ".pdf")), width = 12, height = 8)
  print(
    ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
      geom_bar(stat = "identity", alpha = 0.7) +
      scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
      labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = contrast) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(size = 0.6),
            plot.title = element_text(hjust = 0.5, size = 18),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom") +
      coord_flip()
  )
  dev.off()
  
  # 写入CSV
  write.csv(results_list[[contrast]], file = file.path(csv_dir, paste0(contrast, ".csv")), row.names = TRUE)
}


write.csv(results_list[["NF_vs_NFMac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF_vs_NFMac.csv", row.names = TRUE)
write.csv(results_list[["NF_vs_NF510"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF_vs_NF510.csv", row.names = TRUE)
write.csv(results_list[["NF_vs_NF510Mac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF_vs_NF510Mac.csv", row.names = TRUE)
write.csv(results_list[["NFMac_vs_NF510"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NFMac_vs_NF510.csv", row.names = TRUE)
write.csv(results_list[["NFMac_vs_NF510Mac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NFMac_vs_NF510Mac.csv", row.names = TRUE)
write.csv(results_list[["NF510_vs_NF510Mac"]], file = "~/LB/data/RNA-seq/EN2_510/fig/EN2_510_Mac/GSVA/NF510_vs_NF510Mac.csv", row.names = TRUE)




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

counts <- read.csv(file = '~/LB/data/RNA-seq/BC/co/data/BN_counts_Mac.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型
# 样本名称和对应的条件 (替换 Mac 为 T)
sample <- c("BN1", "BN2", "BN3", "BNMac1", "BNMac2", "BNMac3",
            "BNBCO1", "BNBCO2", "BNBCO3", "BNBCOMac1", "BNBCOMac2", "BNBCOMac3")

# 组标签修改 
group_labels <- c(
  rep("BN", 3),        # 原NF组 
  rep("BNMac", 3),     # 原NFMac组 
  rep("BNBCO", 3),     # 原NF510组 
  rep("BNBCOMac", 3)   # 原NF510Mac组 
)


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
desired_order <- c("BN", "BNMac", "BNBCO", "BNBCOMac")
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
pdf("~/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 14, height = 12)
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


# 定义一个函数来获取每个亚群的前6个通路
get_top6_per_cluster <- function(gsva_results, sample_info) {
  unique_clusters <- unique(sample_info$condition)
  top6_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, sample_info$condition == cluster]
    average_scores <- rowMeans(cluster_expr)
    top6 <- names(sort(average_scores, decreasing = TRUE))[1:6]
    top6_per_cluster[[cluster]] <- top6
  }
  
  return(top6_per_cluster)
}

# 获取每个亚群的前6个通路（所有通路一起考虑）
top6_per_cluster <- get_top6_per_cluster(gsva_results, sample_info)

# 构建新的数据矩阵，包含所有选定的通路
selected_pathways <- unique(unlist(top6_per_cluster))

# 选择对应的 GSVA 评分
selected_scores <- gsva_results[selected_pathways, ]

# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(selected_scores)) {
  stop("The number of samples does not match the number of columns in selected_scores.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("BN", "BNMac", "BNBCO", "BNBCOMac")
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
pdf("~/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/GSVA/GSVA_heatmap_top6_per_cluster.pdf", width = 14, height = 12)
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
  main = "GSVA Scores Heatmap (Top 6 Pathways per Cluster)",
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
output_path <- "~/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)





######BN_BCO_T#####
######两两对比火山图####
library(RColorBrewer)
rdylbu7 <-rev(brewer.pal(7,"RdYlBu"))



library(limma)
library(edgeR)  # for voom function

# 假设counts已经加载
counts <- read.csv(file = '~/LB/data/RNA-seq/BC/co/data/BN_counts_T.csv', header = TRUE, row.names = 1)

# 样本名称修改（NF→BN，510→BCO）
sample <- c("BN1", "BN2", "BN3", "BNT1", "BNT2", "BNT3",
            "BNBCO1", "BNBCO2", "BNBCO3", "BNBCOT1", "BNBCOT2", "BNBCOT3")

# 组标签修改 
group_labels <- c(
  rep("BN", 3),        # 原NF组 
  rep("BNT", 3),     # 原NFMac组 
  rep("BNBCO", 3),     # 原NF510组 
  rep("BNBCOT", 3)   # 原NF510Mac组 
)

# 因子水平修改 
group <- factor(group_labels, levels = c("BN", "BNT", "BNBCO", "BNBCOT"))


design <- model.matrix(~0 + group)
colnames(design) <- levels(group)



counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型

counts[is.na(counts)] <- 0
v <- voom(counts, design, plot = FALSE)

# 拟合线性模型
fit <- lmFit(v, design)

# 对比矩阵修改 
contrasts_matrix <- makeContrasts(
  BN_vs_BNBCO = BNBCO - BN,
  BN_vs_BNT = BNT - BN,
  BN_vs_BNBCOT = BNBCOT - BN,
  BNT_vs_BNBCO = BNBCO - BNT,
  BNBCO_vs_BNBCOT = BNBCOT - BNBCO,
  BNT_vs_BNBCOT = BNBCOT - BNT,
  levels = design 
)
# 应用对比
fit2 <- contrasts.fit(fit, contrasts_matrix)
eb <- eBayes(fit2)

# 检查对比矩阵是否成功应用并且有相应的系数
print(colnames(eb$coefficients))

# 创建一个函数来执行差异表达分析并绘制火山图
analyze_and_plot_limma <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                             box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.5, "lines"),
                             segment.color = "grey50")
  
  return(p)
}

# 对每一对对比执行分析并绘制火山图
for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PDF文件路径
    pdf_path <- file.path("/home/data/t070206/LB/data/RNA-seq/BC/co/fig/BN_BCO_T/Volcano/", paste0("Volcano Plot for ", contrast_name, ".pdf"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    
    # 开始绘制 PDF 文件
    pdf(pdf_path, width = 9, height = 6)
    
    # 打印图形到PDF
    print(plot)
    
    # 关闭PDF设备
    dev.off()
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}

analyze_and_plot_limma1 <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) 
  
  return(p)
}


for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma1(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PNG文件路径
    png_path <- file.path("/home/data/t070206/LB/data/RNA-seq/BC/co/fig/BN_BCO_Mac/Volcano/", paste0("Volcano Plot for ", contrast_name, ".png"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
    
    # 使用 ggsave 保存 PNG 文件，指定宽度和高度以英寸为单位
    ggsave(filename = png_path, plot = plot, width = 9, height = 6, units = "in", dpi = 300)
    
    cat("对比", contrast_name, "已生成并保存为PNG图片\n")
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}




output_dir <- "/home/data/t070206/LB/data/RNA-seq/BC/co/fig/BN_BCO_T/Volcano/"

# 确保输出目录存在，如果不存在则创建
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 对每一对比执行分析并导出为CSV文件
for (contrast_name in colnames(eb$coefficients)) {
  # 获取该对比的差异表达结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 构建完整的文件路径
  file_path <- file.path(output_dir, paste0(contrast_name, "_DEGs.csv"))
  
  # 将结果写入CSV文件
  write.csv(results, file=file_path, row.names=TRUE)
  
  cat("已保存文件:", file_path, "\n")
}



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

counts <- read.csv(file = '~/LB/data/RNA-seq/BC/co/data/BN_counts_T.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型


# 样本名称和对应的条件
sample <- c(
  "BN1", "BN2", "BN3", 
  "BNT1", "BNT2", "BNT3",        # BNMac → BNT 
  "BNBCO1", "BNBCO2", "BNBCO3", 
  "BNBCOT1", "BNBCOT2", "BNBCOT3"  # BNBCOMac → BNBCOT 
)

# 组标签修改 
group_labels <- c(
  rep("BN", 3),        # 基础组（原BN）
  rep("BNT", 3),       # BN + T处理（原BNMac）
  rep("BNBCO", 3),     # BN + BCO处理（原BNBCO）
  rep("BNBCOT", 3)     # BN + BCO + T联合处理（原BNBCOMac）
)

# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
library(edgeR)
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

#4.创建设计矩阵，不含截距项
group <- factor(group_labels, levels = c("BN", "BNT", "BNBCO", "BNBCOT"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵 (替换 Mac 为 T)
contrast.matrix  <- makeContrasts(
  BNvsBNBCO = BNBCO - BN,          # 基础组 vs BCO处理组 
  BNvsBNT = BNT - BN,              # 基础组 vs T处理组（原BNMac → BNT）
  BNvsBNBCOT = BNBCOT - BN,        # 基础组 vs BCO+T联合处理组（原BNBCOMac → BNBCOT）
  BNTvsBNBCO = BNBCO - BNT,        # T处理组 vs BCO处理组 
  BNBCOvsBNBCOT = BNBCOT - BNBCO,  # BCO处理组 vs BCO+T联合处理组 
  BNTvsBNBCOT = BNBCOT - BNT,      # T处理组 vs BCO+T联合处理组 
  levels = design 
)

# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#7. 结果输出
contrast_names <- colnames(contrast.matrix)
contrast_list <- contrast_names  # 用于后续循环
results_list <- lapply(seq_along(contrast_names), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})
names(results_list) <- contrast_names

fig_dir <- "~/LB/data/RNA-seq/BC/co/fig/BN_BCO_T/GSVA/"
csv_dir <- file.path(fig_dir, "GSVA")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 循环处理每组
for (contrast in contrast_list) {
  df <- data.frame(ID = rownames(results_list[[contrast]]), score = results_list[[contrast]]$logFC)
  df <- df[order(df$score), ]
  lowest_10 <- head(df, 10)
  highest_10 <- tail(df, 10)
  top_bottom_df <- rbind(lowest_10, highest_10)  # ✅ 修复这里
  
  # 画图
  pdf(file.path(fig_dir, paste0("GSVA_", contrast, ".pdf")), width = 12, height = 8)
  print(
    ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
      geom_bar(stat = "identity", alpha = 0.7) +
      scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
      labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = contrast) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(size = 0.6),
            plot.title = element_text(hjust = 0.5, size = 18),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom") +
      coord_flip()
  )
  dev.off()
  
  # 写入CSV
  write.csv(results_list[[contrast]], file = file.path(csv_dir, paste0(contrast, ".csv")), row.names = TRUE)
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

counts <- read.csv(file = '~/LB/data/RNA-seq/BC/co/data/BN_counts_T.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型
# 样本名称和对应的条件 (替换 Mac 为 T)
sample <- c("BN1", "BN2", "BN3", "BNT1", "BNT2", "BNT3",
            "BNBCO1", "BNBCO2", "BNBCO3", "BNBCOT1", "BNBCOT2", "BNBCOT3")

# 组标签修改 
group_labels <- c(
  rep("BN", 3),        # 原NF组 
  rep("BNT", 3),     # 原NFMac组 
  rep("BNBCO", 3),     # 原NF510组 
  rep("BNBCOT", 3)   # 原NF510Mac组 
)

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
desired_order <- c("BN", "BNT", "BNBCO", "BNBCOT")
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
pdf("~/LB/data/RNA-seq/BC/co/fig/BN_BCO_T/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 14, height = 12)
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


# 定义一个函数来获取每个亚群的前6个通路
get_top6_per_cluster <- function(gsva_results, sample_info) {
  unique_clusters <- unique(sample_info$condition)
  top6_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, sample_info$condition == cluster]
    average_scores <- rowMeans(cluster_expr)
    top6 <- names(sort(average_scores, decreasing = TRUE))[1:6]
    top6_per_cluster[[cluster]] <- top6
  }
  
  return(top6_per_cluster)
}

# 获取每个亚群的前6个通路（所有通路一起考虑）
top6_per_cluster <- get_top6_per_cluster(gsva_results, sample_info)

# 构建新的数据矩阵，包含所有选定的通路
selected_pathways <- unique(unlist(top6_per_cluster))

# 选择对应的 GSVA 评分
selected_scores <- gsva_results[selected_pathways, ]

# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(selected_scores)) {
  stop("The number of samples does not match the number of columns in selected_scores.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("BN", "BNT", "BNBCO", "BNBCOT")
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
pdf("~/LB/data/RNA-seq/BC/co/fig/BN_BCO_T/GSVA/GSVA_heatmap_top6_per_cluster.pdf", width = 14, height = 12)
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
  main = "GSVA Scores Heatmap (Top 6 Pathways per Cluster)",
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
output_path <- "~/LB/data/RNA-seq/BC/co/fig/BN_BCO_T/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)






######LN_LCO_Mac#####
######两两对比火山图####
library(RColorBrewer)
rdylbu7 <-rev(brewer.pal(7,"RdYlBu"))



library(limma)
library(edgeR)  # for voom function

# 假设counts已经加载
counts <- read.csv(file = '~/LB/data/RNA-seq/LC/co/data/LN_counts_Mac.csv', header = TRUE, row.names = 1)

# 样本名称修改（NF→BN，510→BCO）
sample <- c("LN1", "LN2", "LN3", "LNMac1", "LNMac2", "LNMac3",
            "LNLCO1", "LNLCO2", "LNLCO3", "LNLCOMac1", "LNLCOMac2", "LNLCOMac3")

# 组标签修改 
group_labels <- c(
  rep("LN", 3),        # 原NF组 
  rep("LNMac", 3),     # 原NFMac组 
  rep("LNLCO", 3),     # 原NF510组 
  rep("LNLCOMac", 3)   # 原NF510Mac组 
)

# 因子水平修改 
group <- factor(group_labels, levels = c("LN", "LNMac", "LNLCO", "LNLCOMac"))


design <- model.matrix(~0 + group)
colnames(design) <- levels(group)



counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型
counts[is.na(counts)] <- 0
v <- voom(counts, design, plot = FALSE)

# 拟合线性模型
fit <- lmFit(v, design)

# 对比矩阵修改 
contrasts_matrix <- makeContrasts(
  LN_vs_LNLCO = LNLCO - LN,            # 原BNBCO→LNLCO（BN→LN，BCO→LCO）
  LN_vs_LNMac = LNMac - LN,            # 原BNMac→LNMac（仅替换BN→LN）
  LN_vs_LNLCOMac = LNLCOMac - LN,      # 原BNBCOMac→LNLCOMac（BCO→LCO）
  LNMac_vs_LNLCO = LNLCO - LNMac,      # 原BNMac→LNMac，BNBCO→LNLCO 
  LNLCO_vs_LNLCOMac = LNLCOMac - LNLCO, # 原BNBCOMac→LNLCOMac 
  LNMac_vs_LNLCOMac = LNLCOMac - LNMac, # 原BNBCOMac→LNLCOMac 
  levels = design 
)
# 应用对比
fit2 <- contrasts.fit(fit, contrasts_matrix)
eb <- eBayes(fit2)

# 检查对比矩阵是否成功应用并且有相应的系数
print(colnames(eb$coefficients))

# 创建一个函数来执行差异表达分析并绘制火山图

analyze_and_plot_limma <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                             box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.5, "lines"),
                             segment.color = "grey50")
  
  return(p)
}

# 对每一对对比执行分析并绘制火山图
for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PDF文件路径
    pdf_path <- file.path("/home/data/t070206/LB/data/RNA-seq/LC/co/fig/LN_LCO_Mac/Volcano/", paste0("Volcano Plot for ", contrast_name, ".pdf"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    
    # 开始绘制 PDF 文件
    pdf(pdf_path, width = 9, height = 6)
    
    # 打印图形到PDF
    print(plot)
    
    # 关闭PDF设备
    dev.off()
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}


analyze_and_plot_limma1 <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) 
  
  return(p)
}


for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma1(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PNG文件路径
    png_path <- file.path("/home/data/t070206/LB/data/RNA-seq/LC/co/fig/LN_LCO_Mac/Volcano/", paste0("Volcano Plot for ", contrast_name, ".png"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
    
    # 使用 ggsave 保存 PNG 文件，指定宽度和高度以英寸为单位
    ggsave(filename = png_path, plot = plot, width = 9, height = 6, units = "in", dpi = 300)
    
    cat("对比", contrast_name, "已生成并保存为PNG图片\n")
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}



output_dir <- "/home/data/t070206/LB/data/RNA-seq/LC/co/fig/LN_LCO_Mac/Volcano/"

# 确保输出目录存在，如果不存在则创建
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 对每一对比执行分析并导出为CSV文件
for (contrast_name in colnames(eb$coefficients)) {
  # 获取该对比的差异表达结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 构建完整的文件路径
  file_path <- file.path(output_dir, paste0(contrast_name, "_DEGs.csv"))
  
  # 将结果写入CSV文件
  write.csv(results, file=file_path, row.names=TRUE)
  
  cat("已保存文件:", file_path, "\n")
}






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

counts <- read.csv(file = '~/LB/data/RNA-seq/LC/co/data/LN_counts_Mac.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型


# 样本名称和对应的条件
sample <- c("LN1", "LN2", "LN3", "LNMac1", "LNMac2", "LNMac3",
            "LNLCO1", "LNLCO2", "LNLCO3", "LNLCOMac1", "LNLCOMac2", "LNLCOMac3")

# 组标签修改 
group_labels <- c(
  rep("LN", 3),        # 原NF组 
  rep("LNMac", 3),     # 原NFMac组 
  rep("LNLCO", 3),     # 原NF510组 
  rep("LNLCOMac", 3)   # 原NF510Mac组 
)




# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
library(edgeR)
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

#4.创建设计矩阵，不含截距项
group <- factor(group_labels, levels = c("LN", "LNMac", "LNLCO", "LNLCOMac"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵 (替换 Mac 为 T)
contrast.matrix <- makeContrasts(
  LN_vs_LNLCO = LNLCO - LN,            # 原BNBCO→LNLCO（BN→LN，BCO→LCO）
  LN_vs_LNMac = LNMac - LN,            # 原BNMac→LNMac（仅替换BN→LN）
  LN_vs_LNLCOMac = LNLCOMac - LN,      # 原BNBCOMac→LNLCOMac（BCO→LCO）
  LNMac_vs_LNLCO = LNLCO - LNMac,      # 原BNMac→LNMac，BNBCO→LNLCO 
  LNLCO_vs_LNLCOMac = LNLCOMac - LNLCO, # 原BNBCOMac→LNLCOMac 
  LNMac_vs_LNLCOMac = LNLCOMac - LNMac, # 原BNBCOMac→LNLCOMac 
  levels = design 
)


# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#7. 结果输出
contrast_names <- colnames(contrast.matrix)
contrast_list <- contrast_names  # 用于后续循环
results_list <- lapply(seq_along(contrast_names), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})
names(results_list) <- contrast_names

fig_dir <- "~/LB/data/RNA-seq/LC/co/fig/LN_LCO_Mac/GSVA/"
csv_dir <- file.path(fig_dir, "GSVA")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 循环处理每组
for (contrast in contrast_list) {
  df <- data.frame(ID = rownames(results_list[[contrast]]), score = results_list[[contrast]]$logFC)
  df <- df[order(df$score), ]
  lowest_10 <- head(df, 10)
  highest_10 <- tail(df, 10)
  top_bottom_df <- rbind(lowest_10, highest_10)  # ✅ 修复这里
  
  # 画图
  pdf(file.path(fig_dir, paste0("GSVA_", contrast, ".pdf")), width = 12, height = 8)
  print(
    ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
      geom_bar(stat = "identity", alpha = 0.7) +
      scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
      labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = contrast) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(size = 0.6),
            plot.title = element_text(hjust = 0.5, size = 18),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom") +
      coord_flip()
  )
  dev.off()
  
  # 写入CSV
  write.csv(results_list[[contrast]], file = file.path(csv_dir, paste0(contrast, ".csv")), row.names = TRUE)
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

counts <- read.csv(file = '~/LB/data/RNA-seq/LC/co/data/LN_counts_Mac.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型
# 样本名称和对应的条件 (替换 Mac 为 T)
sample <- c("LN1", "LN2", "LN3", "LNMac1", "LNMac2", "LNMac3",
            "LNLCO1", "LNLCO2", "LNLCO3", "LNLCOMac1", "LNLCOMac2", "LNLCOMac3")

# 组标签修改 
group_labels <- c(
  rep("LN", 3),        # 原NF组 
  rep("LNMac", 3),     # 原NFMac组 
  rep("LNLCO", 3),     # 原NF510组 
  rep("LNLCOMac", 3)   # 原NF510Mac组 
)



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
desired_order <- c("LN", "LNMac", "LNLCO", "LNLCOMac")
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
pdf("~/LB/data/RNA-seq/LC/co/fig/LN_LCO_Mac/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 14, height = 12)
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


# 定义一个函数来获取每个亚群的前6个通路
get_top6_per_cluster <- function(gsva_results, sample_info) {
  unique_clusters <- unique(sample_info$condition)
  top6_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, sample_info$condition == cluster]
    average_scores <- rowMeans(cluster_expr)
    top6 <- names(sort(average_scores, decreasing = TRUE))[1:6]
    top6_per_cluster[[cluster]] <- top6
  }
  
  return(top6_per_cluster)
}

# 获取每个亚群的前6个通路（所有通路一起考虑）
top6_per_cluster <- get_top6_per_cluster(gsva_results, sample_info)

# 构建新的数据矩阵，包含所有选定的通路
selected_pathways <- unique(unlist(top6_per_cluster))

# 选择对应的 GSVA 评分
selected_scores <- gsva_results[selected_pathways, ]

# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(selected_scores)) {
  stop("The number of samples does not match the number of columns in selected_scores.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("LN", "LNMac", "LNLCO", "LNLCOMac")
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
pdf("~/LB/data/RNA-seq/LC/co/fig/LN_LCO_Mac/GSVA/GSVA_heatmap_top6_per_cluster.pdf", width = 14, height = 12)
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
  main = "GSVA Scores Heatmap (Top 6 Pathways per Cluster)",
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
output_path <- "~/LB/data/RNA-seq/LC/co/fig/LN_LCO_Mac/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)





######LN_LCO_T#####
######两两对比火山图####
library(RColorBrewer)
rdylbu7 <-rev(brewer.pal(7,"RdYlBu"))



library(limma)
library(edgeR)  # for voom function

# 假设counts已经加载
counts <- read.csv(file = '~/LB/data/RNA-seq/LC/co/data/LN_counts_T.csv', header = TRUE, row.names = 1)

# 样本名称修改（NF→BN，510→BCO）
sample <- c("LN1", "LN2", "LN3", "LNT1", "LNT2", "LNT3", 
            "LNLCO1", "LNLCO2", "LNLCO3", "LNLCOT1", "LNLCOT2", "LNLCOT3")
group_labels <- c(rep("LN", 3), rep("LNT", 3), rep("LNLCO", 3), rep("LNLCOT", 3))

# 因子水平修改 
group <- factor(group_labels, levels = c("LN", "LNT", "LNLCO", "LNLCOT"))


design <- model.matrix(~0 + group)
colnames(design) <- levels(group)



counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型

counts[is.na(counts)] <- 0
v <- voom(counts, design, plot = FALSE)

# 拟合线性模型
fit <- lmFit(v, design)

# 对比矩阵修改 
contrasts_matrix <- makeContrasts(
  LN_vs_LNLCO = LNLCO - LN,            # 原BNBCO→LNLCO（BN→LN，BCO→LCO）
  LN_vs_LNT = LNT - LN,                # 原BNT→LNT（BN→LN，T保留）
  LN_vs_LNLCOT = LNLCOT - LN,          # 原BNBCOT→LNLCOT（BN→LN，BCO→LCO）
  LNT_vs_LNLCO = LNLCO - LNT,          # 原BNT→LNT，BNBCO→LNLCO 
  LNLCO_vs_LNLCOT = LNLCOT - LNLCO,    # 原BNBCOT→LNLCOT 
  LNT_vs_LNLCOT = LNLCOT - LNT,        # 原BNBCOT→LNLCOT 
  levels = design 
)
# 应用对比
fit2 <- contrasts.fit(fit, contrasts_matrix)
eb <- eBayes(fit2)

# 检查对比矩阵是否成功应用并且有相应的系数
print(colnames(eb$coefficients))

# 创建一个函数来执行差异表达分析并绘制火山图
analyze_and_plot_limma <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene), 
                             box.padding = unit(0.35, "lines"), 
                             point.padding = unit(0.5, "lines"),
                             segment.color = "grey50")
  
  return(p)
}

# 对每一对对比执行分析并绘制火山图
for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PDF文件路径
    pdf_path <- file.path("/home/data/t070206/LB/data/RNA-seq/LC/co/fig/LN_LCO_T/Volcano/", paste0("Volcano Plot for ", contrast_name, ".pdf"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    
    # 开始绘制 PDF 文件
    pdf(pdf_path, width = 9, height = 6)
    
    # 打印图形到PDF
    print(plot)
    
    # 关闭PDF设备
    dev.off()
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}

analyze_and_plot_limma1 <- function(contrast_name, eb) {
  if (!(contrast_name %in% colnames(eb$coefficients))) {
    stop(paste("对比名称", contrast_name, "不在eb对象的系数中"))
  }
  
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  if (is.null(results) || nrow(results) == 0) {
    warning(paste("对比", contrast_name, "没有返回任何结果"))
    return(NULL)
  }
  
  df <- data.frame(
    logFC = results$logFC,
    pval = results$P.Value,
    gene = rownames(results),
    fdr = results$adj.P.Val
  )
  
  df$color <- "grey"
  df$color[df$logFC > 1] <- "#E41A1C"
  df$color[df$logFC < -1] <- "#377EB8"
  
  top_up <- head(df[order(-df$logFC), ], 30)
  top_down <- head(df[order(df$logFC), ], 30)
  top_genes <- rbind(top_up, top_down)
  
  p <- ggplot(df, aes(x = logFC, y = -log10(pval), color = color)) +
    geom_point(alpha = 0.6) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "Log2 Fold Change", y = "-log10(P-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) 
  
  return(p)
}


for (contrast_name in colnames(eb$coefficients)) {  # 使用colnames而不是rownames
  plot <- analyze_and_plot_limma1(contrast_name, eb)
  
  if (!is.null(plot)) {
    # 定义PNG文件路径
    png_path <- file.path("/home/data/t070206/LB/data/RNA-seq/LC/co/fig/LN_LCO_T/Volcano/", paste0("Volcano Plot for ", contrast_name, ".png"))
    
    # 确保输出目录存在，如果不存在则创建
    dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
    
    # 使用 ggsave 保存 PNG 文件，指定宽度和高度以英寸为单位
    ggsave(filename = png_path, plot = plot, width = 9, height = 6, units = "in", dpi = 300)
    
    cat("对比", contrast_name, "已生成并保存为PNG图片\n")
  } else {
    cat("对比", contrast_name, "未生成图表\n")
  }
}




output_dir <- "/home/data/t070206/LB/data/RNA-seq/LC/co/fig/LN_LCO_T/Volcano/"

# 确保输出目录存在，如果不存在则创建
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 对每一对比执行分析并导出为CSV文件
for (contrast_name in colnames(eb$coefficients)) {
  # 获取该对比的差异表达结果
  results <- topTable(eb, coef=contrast_name, number=Inf, adjust="fdr")
  
  # 构建完整的文件路径
  file_path <- file.path(output_dir, paste0(contrast_name, "_DEGs.csv"))
  
  # 将结果写入CSV文件
  write.csv(results, file=file_path, row.names=TRUE)
  
  cat("已保存文件:", file_path, "\n")
}



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

counts <- read.csv(file = '~/LB/data/RNA-seq/LC/co/data/LN_counts_T.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型


# 样本名称和对应的条件
sample <- c("LN1", "LN2", "LN3", "LNT1", "LNT2", "LNT3", 
            "LNLCO1", "LNLCO2", "LNLCO3", "LNLCOT1", "LNLCOT2", "LNLCOT3")
group_labels <- c(rep("LN", 3), rep("LNT", 3), rep("LNLCO", 3), rep("LNLCOT", 3))




# 创建样本信息表
sample_info <- data.frame(sample = sample, condition = group_labels, stringsAsFactors = FALSE)

# 1. 数据准备 - 标准化 counts 数据
# 假设 'counts' 是你的原始计数矩阵 (基因 x 样本)
library(edgeR)
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

#4.创建设计矩阵，不含截距项
group <- factor(group_labels, levels = c("LN", "LNT", "LNLCO", "LNLCOT"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
# 使用voom转换计数数据并创建limma对象
v <- voom(counts, design, plot=FALSE)

# 5. 拟合线性模型
fit <- lmFit(gsva_results, design)

# 6. 定义并应用两两对比矩阵 (替换 Mac 为 T)

contrast.matrix <- makeContrasts(
  LN_vs_LNLCO = LNLCO - LN,            # 原BNBCO→LNLCO（BN→LN，BCO→LCO）
  LN_vs_LNT = LNT - LN,                # 原BNT→LNT（BN→LN，T保留）
  LN_vs_LNLCOT = LNLCOT - LN,          # 原BNBCOT→LNLCOT（BN→LN，BCO→LCO）
  LNT_vs_LNLCO = LNLCO - LNT,          # 原BNT→LNT，BNBCO→LNLCO 
  LNLCO_vs_LNLCOT = LNLCOT - LNLCO,    # 原BNBCOT→LNLCOT 
  LNT_vs_LNLCOT = LNLCOT - LNT,        # 原BNBCOT→LNLCOT 
  levels = design 
)


# 应用对比并估计平均差异
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#7. 结果输出
contrast_names <- colnames(contrast.matrix)
contrast_list <- contrast_names  # 用于后续循环
results_list <- lapply(seq_along(contrast_names), function(i) {
  topTable(fit2, coef = i, number = Inf, adjust = "fdr")
})
names(results_list) <- contrast_names

fig_dir <- "~/LB/data/RNA-seq/LC/co/fig/LN_LCO_T/GSVA/"
csv_dir <- file.path(fig_dir, "GSVA")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# 循环处理每组
for (contrast in contrast_list) {
  df <- data.frame(ID = rownames(results_list[[contrast]]), score = results_list[[contrast]]$logFC)
  df <- df[order(df$score), ]
  lowest_10 <- head(df, 10)
  highest_10 <- tail(df, 10)
  top_bottom_df <- rbind(lowest_10, highest_10)  # ✅ 修复这里
  
  # 画图
  pdf(file.path(fig_dir, paste0("GSVA_", contrast, ".pdf")), width = 12, height = 8)
  print(
    ggplot(top_bottom_df, aes(x = reorder(ID, score), y = score, fill = ifelse(score > 0, "up", "down"))) +
      geom_bar(stat = "identity", alpha = 0.7) +
      scale_fill_manual(values = c("up" = "#E41A1C", "down" = "#377EB8"), name = "Regulation") +
      labs(x = "Pathways", y = "Log Fold Change of GSVA score", title = contrast) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            panel.border = element_rect(size = 0.6),
            plot.title = element_text(hjust = 0.5, size = 18),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom") +
      coord_flip()
  )
  dev.off()
  
  # 写入CSV
  write.csv(results_list[[contrast]], file = file.path(csv_dir, paste0(contrast, ".csv")), row.names = TRUE)
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

counts <- read.csv(file = '~/LB/data/RNA-seq/LC/co/data/LN_counts_T.csv', header = TRUE, row.names = 1)
counts[counts < 0] <- 0 
counts[is.na(counts)] <- 0
# 使用voom转换计数数据并创建limma对象
counts <- as.matrix(counts)  # 转成矩阵
storage.mode(counts) <- "numeric"  # 确保是数值型
# 样本名称修改（NF→BN，510→BCO）
sample <- c("LN1", "LN2", "LN3", "LNT1", "LNT2", "LNT3", 
            "LNLCO1", "LNLCO2", "LNLCO3", "LNLCOT1", "LNLCOT2", "LNLCOT3")
group_labels <- c(rep("LN", 3), rep("LNT", 3), rep("LNLCO", 3), rep("LNLCOT", 3))



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
desired_order <- c("LN", "LNT", "LNLCO", "LNLCOT")
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
pdf("~/LB/data/RNA-seq/LC/co/fig/LN_LCO_T/GSVA/GSVA_heatmap_top12_per_cluster.pdf", width = 14, height = 12)
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


# 定义一个函数来获取每个亚群的前6个通路
get_top6_per_cluster <- function(gsva_results, sample_info) {
  unique_clusters <- unique(sample_info$condition)
  top6_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, sample_info$condition == cluster]
    average_scores <- rowMeans(cluster_expr)
    top6 <- names(sort(average_scores, decreasing = TRUE))[1:6]
    top6_per_cluster[[cluster]] <- top6
  }
  
  return(top6_per_cluster)
}

# 获取每个亚群的前6个通路（所有通路一起考虑）
top6_per_cluster <- get_top6_per_cluster(gsva_results, sample_info)

# 构建新的数据矩阵，包含所有选定的通路
selected_pathways <- unique(unlist(top6_per_cluster))

# 选择对应的 GSVA 评分
selected_scores <- gsva_results[selected_pathways, ]

# 确保 sample_info 的长度与 selected_scores 的列数相匹配
if (nrow(sample_info) != ncol(selected_scores)) {
  stop("The number of samples does not match the number of columns in selected_scores.")
}

# 按照指定顺序重新排列 sample_info
desired_order <- c("LN", "LNT", "LNLCO", "LNLCOT")
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
pdf("~/LB/data/RNA-seq/LC/co/fig/LN_LCO_T/GSVA/GSVA_heatmap_top6_per_cluster.pdf", width = 14, height = 12)
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
  main = "GSVA Scores Heatmap (Top 6 Pathways per Cluster)",
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
output_path <- "~/LB/data/RNA-seq/LC/co/fig/LN_LCO_T/GSVA/GSVA_heatmap_top12_per_cluster.csv"
write.csv(gsva_averages, file = output_path, row.names = FALSE)























