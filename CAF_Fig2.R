

#加载需要用的包
library (Seurat)
library (dplyr)
library (patchwork)
library (sctransform)
library (tidyverse)
library (harmony)
library (future)
library(BPCells)
library(ComplexHeatmap)
library(plot1cell)
library(qs)

#处理iCAF数据
iCAF_rpca_GEO <- readRDS(file = '~/LB/data/pan cancer/CAF/rds/iCAF_rpca_GEO.rds')
Idents(iCAF_rpca_GEO) = "RNA_rpca_snn_res.1.2"
markers = FindAllMarkers(iCAF_rpca_GEO, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.7)
write.csv(markers,file = '~/LB/data/pan cancer/CAF/marker/iCAF_rpca_GEO_1.2.csv')
iCAF_rpca_GEO1 <- subset (iCAF_rpca_GEO,  ident = c("17","22"), invert = TRUE)
Idents(iCAF_rpca_GEO1) <- iCAF_rpca_GEO1$rename
iCAF_rpca_GEO1 <- RenameIdents(iCAF_rpca_GEO1, `Teff_iCAF` = 'apCAF')
iCAF_rpca_GEO1$rename <-Idents(iCAF_rpca_GEO1)
iCAF_rpca_GEO1$rename = factor(iCAF_rpca_GEO1$rename, levels = c("MonoMac_iCAF", "Mac1_iCAF", "Mac2_iCAF", "Gran_iCAF", "TLS_iCAF", "Tem_iCAF", "B_iCAF", "apCAF"))
Idents(iCAF_rpca_GEO1) <- iCAF_rpca_GEO1$rename
saveRDS (iCAF_rpca_GEO1, file = '~/LB/data/pan cancer/CAF/rds/iCAF_rpca_GEO_final.rds')

#Fig 2D
g<-iCAF_rpca_GEO1
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/umap_rpca_GEO.pdf", width = 11, height = 10)
DimPlot(g, label = T, group.by = "rename", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["rename"]][,1]))))) + ggtitle("cluster")
DimPlot(g, label = T, group.by = "RNA_rpca_snn_res.0.1", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_rpca_snn_res.0.1"]][,1]))))) + ggtitle("0.1")
DimPlot(g, label = T, group.by = "RNA_rpca_snn_res.0.3", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_rpca_snn_res.0.3"]][,1]))))) + ggtitle("0.3")
DimPlot(g, label = T, group.by = "RNA_rpca_snn_res.0.4", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_rpca_snn_res.0.4"]][,1]))))) + ggtitle("0.4")
DimPlot(g, label = T, group.by = "RNA_rpca_snn_res.0.6", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_rpca_snn_res.0.6"]][,1]))))) + ggtitle("0.6")
DimPlot(g, label = T, group.by = "RNA_rpca_snn_res.0.8", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_rpca_snn_res.0.8"]][,1]))))) + ggtitle("0.8")
DimPlot(g, label = T, group.by = "RNA_rpca_snn_res.1", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_rpca_snn_res.1"]][,1]))))) + ggtitle("1")
DimPlot(g, label = T, group.by = "RNA_rpca_snn_res.1.2", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = scPalette(length(names(table(g[["RNA_rpca_snn_res.1.2"]][,1]))))) + ggtitle("1.2")
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


umap_data <- Embeddings(g, "umap")


plot_data <- data.frame(
  Global_UMAP_1 = umap_data[, 1],
  Global_UMAP_2 = umap_data[, 2],
  megacluster = g@meta.data$rename,
  UMAP_1 = umap_data[, 1],
  UMAP_2 = umap_data[, 2]
)

ggplot(plot_data %>% arrange(megacluster), aes(x = Global_UMAP_1, y = Global_UMAP_2)) +
  geom_point(aes(color = megacluster), size = .1, alpha = .9) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/umap_rpca_GEO_renameUMAP.pdf", width = 11, height = 10)
ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = megacluster), size = .1, alpha = .8) +
  theme_cowplot(font_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 5, height = 1), ncol = 1))+
  scale_color_manual(values = pal_20)+ggtitle("megacluster")
dev.off()

#count
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/counts cluster.pdf", width = 6, height = 5)
ggplot(g@meta.data, aes(x=`rename`, fill=GEO)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=Tissue)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=Type)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=Gender)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=Stage)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=Pathology)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=Tumor)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=`T`)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=N)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`rename`, fill=M)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/counts cluster_fan.pdf", width = 6, height = 5)
ggplot(g@meta.data, aes(x=`GEO`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`Tissue`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`Type`, fill=rename)) + 
  geom_bar(data = filter(g@meta.data, Type != "NC"), position = "fill") +
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`Gender`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`Stage`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`Pathology`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`Tumor`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`T`, fill=`rename`)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`N`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
ggplot(g@meta.data, aes(x=`M`, fill=rename)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_20) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()
dev.off()

#Fig 2E
markers = FindAllMarkers(g, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.7)
write.csv(markers,file = '~/LB/data/pan cancer/CAF/marker/iCAF_rpca_rename_final.csv')
p1<-FeaturePlot(g, features = c('CCL8','CCL13','CXCL14','CXCR4','CXCL8','CSF2','CSF3','CXCL5','CCL19','CCL21','CXCL9','CXCL10','CD74','HLA-DQA1','IGHA1'))
ggsave("~/LB/data/pan cancer/CAF/fig/iCAF/iCAF_marker.png",plot = p1, width = 16,height = 12,dpi = 300)



#Fig 2F

top10sig <- marker %>%
  group_by(cluster) %>%
  top_n(10, abs(avg_log2FC))


marker$size <- ifelse(marker$gene %in% top10sig$gene, 2, 1)

cluster_colors <- data.frame(
  cluster = c("MonoMac_iCAF", "Mac1_iCAF", "Mac2_iCAF", "Gran_iCAF", "TLS_iCAF", "Tem_iCAF", "B_iCAF", "apCAF"),
  xmin = seq(1, 8, by = 1),
  xmax = seq(1.8, 8.8, by = 1),
  color = pal_20[seq_len(8)]  
)

ymin_val <- min(marker$avg_log2FC) - diff(range(marker$avg_log2FC)) * 0.1 

tile_data <- cluster_colors %>%
  mutate(ymin = ymin_val,
         ymax = ymin_val + diff(range(marker$avg_log2FC)) * 0.1, 
         x = (xmin + xmax) / 2,  
         y = ymin)  


pdf("~/LB/data/pan cancer/CAF/fig/iCAF/多组火山图.pdf", width = 12, height = 6)
ggplot(marker, aes(x = as.numeric(cluster), y = avg_log2FC, color = p_val_adj < 0.01)) + 
  geom_jitter(aes(size = factor(size)), width = 0.4, alpha = 0.8) +
  scale_size_manual(values = c("1" = 0.85, "2" = 1)) + 
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) + 
  geom_text_repel(data = top10sig, aes(label = gene), size = 2.5) + 
  
  # 使用geom_tile来绘制色块
  geom_tile(data = tile_data,
            aes(x = x-0.4, y = 0, fill = color),
            height = diff(range(marker$avg_log2FC)) * 0.1,  # 根据实际数据范围调整高度
            color = "black",
            alpha = 0.6,
            show.legend = FALSE) + 
  
  labs(title = "Volcano Plot with Highlighted Top Genes and Cluster Tiles",
       x = "Cluster",
       y = "log2 Fold Change",
       color = "Significant (adj. P-val)",
       fill = "Cluster") + 
  
  scale_fill_identity(guide = "legend") + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    axis.line.y = element_line(color = "black", size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 15)
  )
dev.off() 

#Fig 2G-J and Fig S8

expr <- GetAssayData(g, slot = "data")
expr <- as.matrix(expr)
rename <- g@meta.data$rename
msgdC2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")

keggSet <- split(msgdC2$gene_symbol, msgdC2$gs_name)

zpar <- zscoreParam(expr, geneSets = keggSet)
gsva_results_kegg <- gsva(zpar, verbose = TRUE)

#gopb
msgdC5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gobpSet <- split(msgdC5$gene_symbol, msgdC5$gs_name)
gobpSet <- lapply(gobpSet, function(set) set[length(set) > 1])
gobpSet <- gobpSet[sapply(gobpSet, length) > 1]
zpar <- zscoreParam(expr, geneSets = gobpSet, minSize = 2)
gsva_results_gobp <- gsva(zpar, verbose = TRUE)

average_scores <- rowMeans(gsva_results_kegg)
top_50_pathways <- names(sort(average_scores, decreasing = TRUE))[1:50]
top_50_scores <- gsva_results_kegg[top_50_pathways, ]
cluster_scores <- aggregate(t(top_50_scores), by = list(cluster = rename), FUN = mean)

heatmap_data <- as.matrix(cluster_scores[, -1])
rownames(heatmap_data) <- cluster_scores$cluster

pdf("~/LB/data/pan cancer/CAF/fig/GSVA/GSVA_KEGG_heatmap.pdf", width = 12, height = 8)
pheatmap(
  heatmap_data,
  scale = "row", 
  cluster_rows = FALSE,  
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap"
)
dev.off()


bmem_expr <- expr[, rename == "Mac2_iCAF"]
glyco_genes <- keggSet[["KEGG_GLYCOLYSIS_GLUCONEOGENESIS"]]
glyco_expr <- bmem_expr[rownames(bmem_expr) %in% glyco_genes, ]
average_glyco_expr <- rowMeans(glyco_expr)

top_glyco_genes <- sort(average_glyco_expr, decreasing = TRUE)
top_genes_glyco <- names(top_glyco_genes)[1:30]
top_gene_expr_glyco <- top_glyco_genes[1:30]

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/Gene_Expression_Mac2_iCAF_Glycolysis-Glucogenesis.pdf", width = 12, height = 8)
barplot(
  top_gene_expr_glyco,
  names.arg = top_genes_glyco,
  col = "steelblue",
  las = 2, 
  main = "Top 30 Highly Expressed Genes in Mac2_iCAF for Glycolysis-Glucogenesis Pathway",
  ylab = "Average Expression"
)
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/GSVA/GSVA_GOBP_heatmap.pdf", width = 12, height = 15)
pheatmap(
  heatmap_data,
  scale = "row", 
  cluster_rows = FALSE, 
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap for GOBP")
dev.off()


msgdH <- msigdbr(species = "Homo sapiens", category = "H")
hallmarkSet <- split(msgdH$gene_symbol, msgdH$gs_name)
zpar <- zscoreParam(expr, geneSets = hallmarkSet)
gsva_results_hallmark <- gsva(zpar, verbose = TRUE)
average_scores <- rowMeans(gsva_results)
top_40_pathways <- names(sort(average_scores, decreasing = TRUE))[1:40]
top_40_scores <- gsva_results[top_40_pathways, ]

cluster_scores <- aggregate(t(top_40_scores), by = list(cluster = rename), FUN = mean)
heatmap_data <- as.matrix(cluster_scores[, -1])
rownames(heatmap_data) <- cluster_scores$cluster
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/GSVA/GSVA_Hallmark_heatmap.pdf", width = 12, height = 15)
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = FALSE,  
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap for Hallmark Pathways"
)
dev.off()



mac2_icaf_expr <- expr[, rename == "Mac2_iCAF"]
mtorc1_genes <- hallmarkSet[["HALLMARK_MTORC1_SIGNALING"]]
mtorc1_expr <- mac2_icaf_expr[rownames(mac2_icaf_expr) %in% mtorc1_genes, ]
average_mtorc1_expr <- rowMeans(mtorc1_expr)
top_mtorc1_genes <- sort(average_mtorc1_expr, decreasing = TRUE)
top_genes_mtorc1 <- names(top_mtorc1_genes)[1:30]
top_gene_expr_mtorc1 <- top_mtorc1_genes[1:30]

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/GSVA/Gene_Expression_Mac2_iCAF_MTORC1.pdf", width = 12, height = 8)
barplot(
  top_gene_expr_mtorc1,
  names.arg = top_genes_mtorc1,
  col = "steelblue",
  las = 2, 
  main = "Top 30 Highly Expressed Genes in Mac2_iCAF for MTORC1 Pathway",
  ylab = "Average Expression"
)
dev.off()
gran_icaf_expr <- expr[, rename == "Gran_iCAF"]

tnfa_nfkB_genes <- hallmarkSet[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]
tnfa_nfkB_expr <- gran_icaf_expr[rownames(gran_icaf_expr) %in% tnfa_nfkB_genes, ]

average_tnfa_nfkB_expr <- rowMeans(tnfa_nfkB_expr)
top_tnfa_nfkB_genes <- sort(average_tnfa_nfkB_expr, decreasing = TRUE)
top_genes_tnfa_nfkB <- names(top_tnfa_nfkB_genes)[1:30]
top_gene_expr_tnfa_nfkB <- top_tnfa_nfkB_genes[1:30]

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/GSVA/Gene_Expression_Gran_iCAF_TNFA-NFKB.pdf", width = 12, height = 8)
barplot(
  top_gene_expr_tnfa_nfkB,
  names.arg = top_genes_tnfa_nfkB,
  col = "steelblue",
  las = 2,  # 横轴标签旋转
  main = "Top 30 Highly Expressed Genes in Gran_iCAF for TNFA-NFKB Pathway",
  ylab = "Average Expression"
)
dev.off()


###整合GSVA结果到seurat对象中##
# 将GSVA结果转置，以便整合到Seurat对象中
gsva_results_kegg_t <- t(gsva_results_kegg)
gsva_results_gobp_t <- t(gsva_results_gobp)
gsva_results_hallmark_t <- t(gsva_results_hallmark)
# 将GSVA结果添加到Seurat对象的元数据中
g <- AddMetaData(g, metadata = as.data.frame(gsva_results_kegg_t))
g <- AddMetaData(g, metadata = as.data.frame(gsva_results_gobp_t))
g <- AddMetaData(g, metadata = as.data.frame(gsva_results_hallmark_t))

p1<-FeaturePlot(g, features = c("KEGG_GLYCOLYSIS_GLUCONEOGENESIS","KEGG_PRIMARY_IMMUNODEFICIENCY","KEGG_COMPLEMENT_AND_COAGULATION_CASCADES","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","GOBP_IMMUNE_RESPONSE","HALLMARK_GLYCOLYSIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INFLAMMATORY_RESPONSE"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
ggsave("~/LB/data/pan cancer/CAF/fig/iCAF/GSVA/iCAF_gsva_pathways_1.png",plot = p1, width = 20,height = 12,dpi = 300)

p2<-FeaturePlot(g, features = c("KEGG_GLYCOLYSIS_GLUCONEOGENESIS","HALLMARK_HYPOXIA","HALLMARK_OXIDATIVE_PHOSPHORYLATION","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION","GOBP_IMMUNE_RESPONSE","HALLMARK_GLYCOLYSIS","HALLMARK_MTORC1_SIGNALING","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INFLAMMATORY_RESPONSE"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
ggsave("~/LB/data/pan cancer/CAF/fig/iCAF/GSVA/iCAF_gsva_pathways_2.png",plot = p2, width = 20,height = 12,dpi = 300)

FeaturePlot(g, features = c("HALLMARK_GLYCOLYSIS"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("KEGG_GLYCOLYSIS_GLUCONEOGENESIS"), reduction = "umap", min.cutoff = 0,  max.cutoff = 10)
FeaturePlot(g, features = c("HALLMARK_MTORC1_SIGNALING"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("HALLMARK_INTERFERON_GAMMA_RESPONSE"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("HALLMARK_INTERFERON_ALPHA_RESPONSE"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("KEGG_PRIMARY_IMMUNODEFICIENCY"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("HALLMARK_INFLAMMATORY_RESPONSE"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"), reduction = "umap", min.cutoff = 0,  max.cutoff = 15)
FeaturePlot(g, features = c("HALLMARK_HYPOXIA"), reduction = "umap", min.cutoff = 0,  max.cutoff = 10)
FeaturePlot(g, features = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION"), reduction = "umap", min.cutoff = 0.5,  max.cutoff = 10)



cell_paths <- list(
  "Mac2_iCAF" = c("KEGG_GLYCOLYSIS_GLUCONEOGENESIS"),
  "apCAF" = c(
    "KEGG_PRIMARY_IMMUNODEFICIENCY",
    "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
    "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
    "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
    "GOBP_IMMUNE_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE"
  ),
  "Tem_iCAF" = c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  )
)


for (cell_type in names(cell_paths)) {

  cell_expr <- expr[, rename == cell_type]
  

  for (pathway in cell_paths[[cell_type]]) {

    if (startsWith(pathway, "KEGG")) {
      pathway_genes <- keggSet[[pathway]]  
    } else if (startsWith(pathway, "GOBP")) {
      pathway_genes <- gobpSet[[pathway]]  
    } else if (startsWith(pathway, "HALLMARK")) {
      pathway_genes <- hallmarkSet[[pathway]]  
    }
    

    pathway_expr <- cell_expr[rownames(cell_expr) %in% pathway_genes, ]
    
  
    average_pathway_expr <- rowMeans(pathway_expr)
    

    top_pathway_genes <- sort(average_pathway_expr, decreasing = TRUE)
    top_genes_pathway <- names(top_pathway_genes)[1:30]
    top_gene_expr_pathway <- top_pathway_genes[1:30]
    

    pdf(paste0("~/LB/data/pan cancer/CAF/fig/", cell_type, "_", pathway, "_Expression.pdf"), width = 12, height = 8)
    barplot(
      top_gene_expr_pathway,
      names.arg = top_genes_pathway,
      col = "steelblue",
      las = 2,  
      main = paste("Top 30 Highly Expressed Genes in", cell_type, "for", pathway),
      ylab = "Average Expression"
    )
    dev.off()
  }
}

saveRDS (g, file = '~/LB/data/pan cancer/CAF/rds/iCAF_rpca_GEO_final_withGSVA.rds')



######Fig 2A and Fig S5#####
Allcell <- readRDS(file = '~/LB/data/pan cancer/first data file/outerBPCells/Allcell.rds')
CAF_harmony_samples_rename <- readRDS(file = '~/LB/data/pan cancer/CAF/rds/CAF_harmony_samples_rename.rds')

Idents(Allcell) <-Allcell$cluster
Allcell@meta.data$cluster1 = Allcell@meta.data$cluster
meta = CAF_harmony_samples_rename@meta.data
meta$rename <- as.character(meta$rename)
Allcell@meta.data$cluster1 <- as.character(Allcell@meta.data$cluster1)
Allcell@meta.data[rownames(meta), "cluster1"] <- meta$rename
Allcell@meta.data$cluster1 = factor(Allcell@meta.data$cluster1, levels = c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW', 'NF', 'sCAF', 'iCAF', 'myCAF', 'vCAF', 'dCAF', 'lsCAF', 'Mural'))
DimPlot(Allcell, reduction = 'full.UMAP_HARM', group.by = 'cluster1')

Allcell[['cluster2']] <- 'none'
Allcell$cluster2[grep('B', Allcell$cluster1)] <- 'B'
Allcell$cluster2[grep('CD4T', Allcell$cluster1)] <- 'CD4T'
Allcell$cluster2[grep('CD8T', Allcell$cluster1)] <- 'CD8T'
Allcell$cluster2[grep('Mac', Allcell$cluster1)] <- 'Mac'
Allcell$cluster2[grep('Mast', Allcell$cluster1)] <- 'Mast'
Allcell$cluster2[grep('mDC', Allcell$cluster1)] <- 'mDC'
Allcell$cluster2[grep('Plas', Allcell$cluster1)] <- 'Plas'
Allcell$cluster2[grep('Endo', Allcell$cluster1)] <- 'Endo'
Allcell$cluster2[grep('Epi', Allcell$cluster1)] <- 'Epi'
Allcell$cluster2[grep('SCW', Allcell$cluster1)] <- 'SCW'

Allcell$cluster2[grep('NF', Allcell$cluster1)] <- 'NF'
Allcell$cluster2[grep('sCAF', Allcell$cluster1)] <- 'sCAF'
Allcell$cluster2[grep('iCAF', Allcell$cluster1)] <- 'iCAF'
Allcell$cluster2[grep('myCAF', Allcell$cluster1)] <- 'myCAF'
Allcell$cluster2[grep('vCAF', Allcell$cluster1)] <- 'vCAF'
Allcell$cluster2[grep('dCAF', Allcell$cluster1)] <- 'dCAF'
Allcell$cluster2[grep('lsCAF', Allcell$cluster1)] <- 'lsCAF'
Allcell$cluster2[grep('Mural', Allcell$cluster1)] <- 'Mural'

Allcell@meta.data$cluster2 = factor(Allcell@meta.data$cluster2, levels = c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW', 'NF', 'sCAF', 'iCAF', 'myCAF', 'vCAF', 'dCAF', 'lsCAF', 'Mural','none'))
Idents(Allcell) <-Allcell$cluster2
g <- subset (Allcell,   ident = c("none"), invert = TRUE)
g = SketchData(object = g, ncells = 200000, method = "LeverageScore", sketched.assay = "sketch")
g[['sketch2']] <- 'All'
g$sketch2[colnames(g@assays$sketch)] <- 'sketch'
g1 <- subset(g, sketch2 == 'sketch')
DefaultAssay(g1) <- "RNA"
g1[['sketch']] <- NULL
g1[["RNA"]]$counts <- as(g1[["RNA"]]$counts, "dgCMatrix")
g1[["RNA"]]$data <- as(g1[["RNA"]]$data, "dgCMatrix")


Idents(g1) = "cluster1"
CellChatDB = CellChatDB.human 
CellChatDB = subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat = createCellChat(object = GetAssayData(g1, layer = 'data', assay = "RNA"), meta = g1@meta.data, group.by = 'cluster1')
cellchat@DB = CellChatDB
cellchat@meta$cluster1 = cellchat@meta$cluster1 %>% droplevels()
cellchat@idents = cellchat@idents %>% droplevels()
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = computeCommunProb(cellchat)
cellchat = filterCommunication(cellchat, min.cells = 10)
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat = computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat = netEmbedding(cellchat, type = "functional", umap.method = "uwot")
cellchat = netClustering(cellchat, type = "functional")
cellchat = computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat = netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat = netClustering(cellchat, type = "structural")

groupSize = as.numeric(table(cellchat@idents))
pdf("~/LB/data/pan cancer/CAF/fig/cellchat/netVisual.circle.pdf", width = 10, height = 10)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat = cellchat@net$weight
pdf("~/LB/data/pan cancer/CAF/fig/cellchat/netVisual.circle.splitted.pdf", width = 20, height = 16)
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathways.show = cellchat@netP$pathway
dir.create("~/LB/data/pan cancer/CAF/fig/cellchat/pathway")
for(i in pathways.show) {
  pdf(paste0("~/LB/data/pan cancer/CAF/fig/cellchat/pathway/", i, ".pdf"), width = 8, height = 8)
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = i, layout = "circle")
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = i, layout = "chord")
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat, signaling = i, color.heatmap = "Reds"))
  dev.off()
}

gg1 = netVisual_heatmap(cellchat, color.heatmap = "Reds")
gg2 = netVisual_heatmap(cellchat, color.heatmap = "Reds",measure = "weight")
pdf("~/LB/data/pan cancer/CAF/fig/cellchat/netVisual.heatmap.pdf", width = 10, height = 5)
gg1 + gg2
dev.off()

dir.create("~/LB/data/pan cancer/CAF/fig/cellchat/LR")
pathways.show.all = cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  gg = netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("~/LB/data/pan cancer/CAF/fig/cellchat/LR/", pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 2, units = 'in', dpi = 300)
}

pdf("~/LB/data/pan cancer/CAF/fig/cellchat/Chord diagram DE_LRP.pdf", width = 15, height = 15)
print(netVisual_aggregate(cellchat, , signaling = pathways.show, layout = "chord", title = "All"))
for (i in pathways.show) {
  print(netVisual_aggregate(cellchat, signaling = i, , layout = "chord", title = i))
}
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/cellchat/all.LRP.bubbleplot.pdf", width = 15, height = 15)
netVisual_bubble(cellchat, sources.use = c(1:18), targets.use = c(1:18), remove.isolate = F)
dev.off()

cell_clusters <-levels(cellchat@idents)
pdf("~/LB/data/pan cancer/CAF/fig/cellchat/all.in.LRP.bubbleplot.pdf", width = 10, height = 10)
for (i in cell_clusters) {
  print(netVisual_bubble(cellchat, sources.use = c(1:18), targets.use = i, remove.isolate = F, title = i))
}
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/cellchat/all.out.LRP.bubbleplot.pdf", width = 10, height = 10)
for (i in cell_clusters) {
  print(netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:18), remove.isolate = F, title = i))
}
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/cellchat/all.LRP.chordplot.pdf", width = 15, height = 15)
netVisual_chord_gene(cellchat, sources.use = c(11:18), targets.use = c(1:10), lab.cex = 0.5, legend.pos.y = 30)
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/cellchat/netAnalysis_signalingRole_network.pdf", width = 8, height = 3)
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.2, font.size = 10)
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/cellchat/netAnalysis_signalingRole_scatter.pdf", width = 6, height = 6)
print(netAnalysis_signalingRole_scatter(cellchat, title = "All"))
for (i in pathways.show) {
  print(netAnalysis_signalingRole_scatter(cellchat, signaling = i, title = i))
}
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 8, height = 10)
ht3 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", width = 8, height = 10)
pdf("~/LB/data/pan cancer/CAF/fig/cellchat/netAnalysis_signalingRole_heatmap.pdf", width = 15, height = 6)
ht1 + ht2 + ht3
dev.off()


pdf("~/LB/data/pan cancer/CAF/fig/cellchat/rank signal pathway barplot.pdf", width = 5, height = 4)
rankNet(cellchat, mode = "single", stacked = T, do.stat = TRUE)
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/cellchat/rank signal pathway barplot per CAF cluster.pdf", width = 5, height = 4)
for (i in c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW', 'NF', 'sCAF', 'iCAF', 'myCAF', 'vCAF', 'dCAF', 'lsCAF', 'Mural')) {
  print(rankNet(cellchat, mode = "single", stacked = T, do.stat = TRUE, sources.use = i) + ggtitle(i))
}
dev.off()




pdf("~/LB/data/pan cancer/CAF/fig/cellchat/Chord diagram DE_LRP CAF sender receiver by DEG method.pdf", width = 8, height = 8)
for (i in c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW', 'NF', 'sCAF', 'iCAF', 'myCAF', 'vCAF', 'dCAF', 'lsCAF', 'Mural')) 
{print(netVisual_chord_gene(cellchat, sources.use = i, targets.use = c(1:10), lab.cex = 0.5,legend.pos.y = 30,title.name = paste0("Source: ", i)) }
dev.off()
  
saveRDS(cellchat, "~/LB/data/pan cancer/CAF/rds/CAF_cellchat.rds")



######Fig 2B-C and Fig S6#####

iCAF_rpca_GEO_final <- readRDS(file = '~/LB/data/pan cancer/CAF/rds/iCAF_rpca_GEO_final.rds')
saveRDS(iCAF_rpca_GEO_final,file = '~/LB/data/pan cancer/CAF/rds/iCAF_rpca_GEO_final.rds')
Allcell <- readRDS(file = '~/LB/data/pan cancer/first data file/outerBPCells/Allcell.rds')
DefaultAssay(Allcell) <- "RNA"
Allcell[['sketch']] <- NULL

Allcell@meta.data$cluster = plyr::mapvalues(Allcell@meta.data$cluster_full, from = c(0:20), to = c("CAF","Epi","Mac","B","CD4T","Epi","CD8T","Epi","Endo","Epi","mDC","Mural","cyc","Epi","Plas","Epi","Mast","Epi","SCW","Epi","Epi"))
Idents(Allcell) <-Allcell$cluster
DimPlot(Allcell, reduction = 'full.UMAP_HARM', group.by = 'cluster')


Allcell@meta.data$cluster3 = Allcell@meta.data$cluster
meta = iCAF_rpca_GEO_final@meta.data
meta$rename <- as.character(meta$rename)
Allcell@meta.data$cluster3 <- as.character(Allcell@meta.data$cluster3)
Allcell@meta.data[rownames(meta), "cluster3"] <- meta$rename
Allcell@meta.data$cluster3 = factor(Allcell@meta.data$cluster3, levels = c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW', 'Mono_iCAF', 'Mac1_iCAF', 'Mac2_iCAF', 'Gran_iCAF', 'TLS_iCAF', 'Tem_iCAF', 'B_iCAF', 'apCAF'))


Allcell[['cluster4']]<- 'none'
Allcell$cluster4[grep('B', Allcell$cluster3)] <- 'B'
Allcell$cluster4[grep('CD4T', Allcell$cluster3)] <- 'CD4T'
Allcell$cluster4[grep('CD8T', Allcell$cluster3)] <- 'CD8T'
Allcell$cluster4[grep('Mac', Allcell$cluster3)] <- 'Mac'
Allcell$cluster4[grep('Mast', Allcell$cluster3)] <- 'Mast'
Allcell$cluster4[grep('mDC', Allcell$cluster3)] <- 'mDC'
Allcell$cluster4[grep('Plas', Allcell$cluster3)] <- 'Plas'
Allcell$cluster4[grep('Endo', Allcell$cluster3)] <- 'Endo'
Allcell$cluster4[grep('Epi', Allcell$cluster3)] <- 'Epi'
Allcell$cluster4[grep('SCW', Allcell$cluster3)] <- 'SCW'
Allcell$cluster4[grep('Mono_iCAF', Allcell$cluster3)] <- 'Mono_iCAF'
Allcell$cluster4[grep('Mac1_iCAF', Allcell$cluster3)] <- 'Mac1_iCAF'
Allcell$cluster4[grep('Mac2_iCAF', Allcell$cluster3)] <- 'Mac2_iCAF'
Allcell$cluster4[grep('Gran_iCAF', Allcell$cluster3)] <- 'Gran_iCAF'
Allcell$cluster4[grep('TLS_iCAF', Allcell$cluster3)] <- 'TLS_iCAF'
Allcell$cluster4[grep('Tem_iCAF', Allcell$cluster3)] <- 'Tem_iCAF'
Allcell$cluster4[grep('B_iCAF', Allcell$cluster3)] <- 'B_iCAF'
Allcell$cluster4[grep('apCAF', Allcell$cluster3)] <- 'apCAF'
Idents(Allcell) <-Allcell$cluster4

g <- subset (Allcell,  cluster4 = c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW', 'Mono_iCAF', 'Mac1_iCAF', 'Mac2_iCAF', 'Gran_iCAF', 'TLS_iCAF', 'Tem_iCAF', 'B_iCAF', 'apCAF'))
g <- subset (Allcell,   ident = c("none"), invert = TRUE)
saveRDS (Allcell, file = '~/LB/data/pan cancer/first data file/outerBPCells/Allcell.rds')

g = SketchData(object = g, ncells = 200000, method = "LeverageScore", sketched.assay = "sketch")
g[['sketch2']] <- 'All'
g$sketch2[colnames(g@assays$sketch)] <- 'sketch'
g1 <- subset(g, sketch2 == 'sketch')
DefaultAssay(g1) <- "RNA"
g1[['sketch']] <- NULL
g1[["RNA"]]$counts <- as(g1[["RNA"]]$counts, "dgCMatrix")
g1[["RNA"]]$data <- as(g1[["RNA"]]$data, "dgCMatrix")


Idents(g1) = "cluster3"
CellChatDB = CellChatDB.human 
CellChatDB = subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat = createCellChat(object = GetAssayData(g1, layer = 'data', assay = "RNA"), meta = g1@meta.data, group.by = 'cluster3')
cellchat@DB = CellChatDB
cellchat@meta$cluster3 = cellchat@meta$cluster3 %>% droplevels()
cellchat@idents = cellchat@idents %>% droplevels()
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = computeCommunProb(cellchat)
cellchat = filterCommunication(cellchat, min.cells = 10)
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat = computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat = netEmbedding(cellchat, type = "functional", umap.method = "uwot")
cellchat = netClustering(cellchat, type = "functional")
cellchat = computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat = netEmbedding(cellchat, type = "structural", umap.method = "uwot")
cellchat = netClustering(cellchat, type = "structural")

groupSize = as.numeric(table(cellchat@idents))
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/netVisual.circle.pdf", width = 10, height = 10)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat = cellchat@net$weight
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/netVisual.circle.splitted.pdf", width = 20, height = 16)
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] = mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathways.show = cellchat@netP$pathway
dir.create("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/pathway")
for(i in pathways.show) {
  pdf(paste0("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/pathway/", i, ".pdf"), width = 8, height = 8)
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = i, layout = "circle")
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = i, layout = "chord")
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat, signaling = i, color.heatmap = "Reds"))
  dev.off()
}

gg1 = netVisual_heatmap(cellchat, color.heatmap = "Reds")
gg2 = netVisual_heatmap(cellchat, color.heatmap = "Reds",measure = "weight")
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/netVisual.heatmap.pdf", width = 10, height = 5)
gg1 + gg2
dev.off()

dir.create("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/LR")
pathways.show.all = cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  gg = netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/LR/", pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 2, units = 'in', dpi = 300)
}

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/Chord diagram DE_LRP.pdf", width = 15, height = 15)
print(netVisual_aggregate(cellchat, , signaling = pathways.show, layout = "chord", title = "All"))
for (i in pathways.show) {
  print(netVisual_aggregate(cellchat, signaling = i, , layout = "chord", title = i))
}
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/all.LRP.bubbleplot.pdf", width = 15, height = 15)
netVisual_bubble(cellchat, sources.use = c(1:18), targets.use = c(1:18), remove.isolate = F)
dev.off()

cell_clusters <-levels(cellchat@idents)
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/all.in.LRP.bubbleplot.pdf", width = 10, height = 10)
for (i in cell_clusters) {
  print(netVisual_bubble(cellchat, sources.use = c(1:18), targets.use = i, remove.isolate = F, title = i))
}
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/all.out.LRP.bubbleplot.pdf", width = 10, height = 10)
for (i in cell_clusters) {
  print(netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:18), remove.isolate = F, title = i))
}
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/all.LRP.chordplot.pdf", width = 15, height = 15)
netVisual_chord_gene(cellchat, sources.use = c(11:18), targets.use = c(1:10), lab.cex = 0.5, legend.pos.y = 30)
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/netAnalysis_signalingRole_network.pdf", width = 8, height = 3)
netAnalysis_signalingRole_network(cellchat, width = 8, height = 2.2, font.size = 10)
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/netAnalysis_signalingRole_scatter.pdf", width = 6, height = 6)
print(netAnalysis_signalingRole_scatter(cellchat, title = "All"))
for (i in pathways.show) {
  print(netAnalysis_signalingRole_scatter(cellchat, signaling = i, title = i))
}
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 8, height = 10)
ht3 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", width = 8, height = 10)
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/netAnalysis_signalingRole_heatmap.pdf", width = 15, height = 6)
ht1 + ht2 + ht3
dev.off()


pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/rank signal pathway barplot.pdf", width = 5, height = 4)
rankNet(cellchat, mode = "single", stacked = T, do.stat = TRUE)
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/rank signal pathway barplot per CAF cluster.pdf", width = 5, height = 4)
for (i in c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW', 'Mono_iCAF', 'Mac1_iCAF', 'Mac2_iCAF', 'Gran_iCAF', 'TLS_iCAF', 'Tem_iCAF', 'B_iCAF', 'apCAF')) {
  print(rankNet(cellchat, mode = "single", stacked = T, do.stat = TRUE, sources.use = i) + ggtitle(i))
}
dev.off()




pdf("~/LB/data/pan cancer/CAF/fig/iCAF/cellchat/Chord diagram DE_LRP CAF sender receiver by DEG method.pdf", width = 8, height = 8)
for (i in c('B', 'CD4T', 'CD8T', 'Mac', 'Mast', 'mDC', 'Plas', 'Endo', 'Epi', 'SCW','Mono_iCAF', 'Mac1_iCAF', 'Mac2_iCAF', 'Gran_iCAF', 'TLS_iCAF', 'Tem_iCAF', 'B_iCAF', 'apCAF')) {
  print(netVisual_chord_gene(cellchat,sources.use = i,targets.use = c(1:10),lab.cex = 0.5, legend.pos.y = 30,title.name = paste0("Source: ", i)))
}
dev.off()







saveRDS(cellchat, "~/LB/data/pan cancer/CAF/rds/iCAF_cellchat.rds")


ave = AverageExpression(g, assays = "RNA", features = VariableFeatures(g, assay = "RNA"), slot = "data", group.by = "cluster")$RNA
ave = as.matrix(ave)
data = cor(ave)
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/corrplot cluster.pdf", width = 6.2, height = 6)
corrplot::corrplot(data, col = c(colorRampPalette(c("navy"))(1000), colorRampPalette(c("navy","white","peachpuff","deeppink4"))(1000)))
dev.off()

rownames_data <- rownames(data)
colnames_data <- colnames(data)

selected_data <- data[rownames_data[(length(rownames_data) - 7):length(rownames_data)], 
                      colnames_data[1:10]]
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/corrplot cluster1.pdf", width = 6.2, height = 6)
corrplot::corrplot(selected_data, col = c( colorRampPalette(c("#00008B", "white","peachpuff", "#FF0000"))(1000)),cl.lim = c(-0.5, 0.5), tl.col = "black", tl.srt = 45, diag = FALSE)
dev.off()



pdf("~/LB/data/pan cancer/CAF/fig/iCAF/ggcorrplot cluster.pdf", width = 7, height = 6)
ggcorrplot::ggcorrplot(data, colors = c("#00008B", "white", "#FF0000"))
dev.off()
pdf("~/LB/data/pan cancer/CAF/fig/iCAF/ggcorrplot cluster1.pdf", width = 7, height = 6)
ggcorrplot::ggcorrplot(selected_data, colors = c("#00008B", "white", "#FF0000"))
dev.off()

pheatmap::pheatmap(selected_data,cluster_rows = FALSE,cluster_cols = FALSE)

pdf("~/LB/data/pan cancer/CAF/fig/iCAF/pheatmap cluster.pdf", width = 12, height = 10)
pheatmap(
  selected_data,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = FALSE,  
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "black",cellwidth = 50, cellheight = 50,
  display_numbers = TRUE,        
  number_color = "black",         
  fontsize=10,                    
  number_format = "%.1f", 
  main = "Corrplot"
)
dev.off()
