#####CAF补充代码

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


##Fig 1B##
library(ggplot2)

tumor_data <- data.frame(
  Tumor_Type = c("HNSCC", "CSCC", "ESCC", "LUAC", "LUSC", "BRCA", "STAD", 
                 "HCC", "ICCA", "CHC", "PDAC", "CRC"),
  Sample_Count = c(18, 10, 60, 10, 5, 26, 31, 40, 12, 20, 35, 63)
)

pdf("~/LB/data/pan cancer/CAF/fig/Fig1/Patients Count.pdf", width = 12, height = 6)
ggplot(tumor_data, aes(x = Tumor_Type, y = Sample_Count)) +
  geom_bar(stat = "identity", fill = "#377EB8", color = "#377EB8", width = 0.8) + 
  theme_minimal() + 
  xlab("Tumor Type") + 
  ylab("Number of Patients") +
  ggtitle("Patients Counts by Tumor Type") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black", size = 1), 
    axis.ticks = element_line(color = "black"), 
    axis.text.y = element_text(size = 12), 
    plot.title = element_text(size = 14, face = "bold"), 
    plot.margin = margin(t = 20, r = 10, b = 30, l = 10, unit = "pt") 
  ) +
  geom_text(aes(label = Sample_Count), vjust = -0.5, color = "black", size = 5) 

dev.off() 

tumor_data <- data.frame(
  Tumor_Type = c("HNSCC", "CSCC", "ESCC", "LUAC", "LUSC", "BRCA", "STAD", 
                 "HCC", "ICCA", "CHC", "PDAC", "CRC"),
  Sample_Count = c(18, 20, 64, 10, 5, 26, 48, 45, 16, 23, 35, 189)
)


p <- ggplot(tumor_data, aes(x = Tumor_Type, y = Sample_Count)) +
  geom_bar(stat = "identity", fill = "#377EB8", color = "#377EB8", width = 0.8) +
  theme_minimal() + 
  xlab("Tumor Type") + 
  ylab("Number of samples") + 
  ggtitle("Sample Counts by Tumor Type") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 20, r = 10, b = 30, l = 10, unit = "pt") 
  ) +
  geom_text(aes(label = Sample_Count), vjust = -0.5, color = "black", size = 5) 


pdf("~/LB/data/pan cancer/CAF/fig/Fig1/Sample Count.pdf", width = 12, height = 6)
print(p)
dev.off()

##Fig S1B##
#CAF Ratio
caf_data <- data.frame(
  Tumor_Type = c("HNSCC", "CSCC", "ESCC", "NSCLC1", "NSCLC2", "BC", "GC", "PLC1", "PLC2", "PDAC", "CRC"),
  CAF_Ratio = c(3.94, 2, 13.36, 2.26, 2.77, 9.19, 8.16, 4.37, 3.86, 11.91, 12.52),
  Source = c("NC_2021", "Cell_2020", "NC_2021", "Immunity_2019", "NM_2018", "NG_2021", "CD_2022", "Cell_2020", "JHep_2021", "CR_2019", "NG_2022")
)

pdf("~/LB/data/pan cancer/CAF/fig/Fig1/CAF Ratio.pdf", width = 12, height = 6)
ggplot(caf_data, aes(x = Tumor_Type, y = CAF_Ratio)) +
  geom_bar(stat = "identity", fill = "#377EB8", color = "#377EB8", width = 0.8) + 
  theme_minimal() + 
  xlab("Tumor Type") + 
  ylab("CAF Ratio (%)") + 
  ggtitle("CAF Ratios by Tumor Type") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 20, r = 10, b = 30, l = 10, unit = "pt")
  ) +
  geom_text(aes(label = sprintf("%.2f", CAF_Ratio)), vjust = -0.5, color = "black", size = 5)

dev.off()

#Number of Cells
cell_subtypes <- data.frame(
  Cell_Type = c("CAF", "Epi", "Mac", "B", "CD4T", "CD8T", "Endo", "mDC", "Mural", "cyc", "Plas", "Mast", "SCW"),
  Cell_Count = c(102297, 212297, 132908, 72271, 231597, 199952, 65811, 67369, 31382, 34199, 11987, 9641, 5205)
)


pdf("~/LB/data/pan cancer/CAF/fig/Fig1/Cell Subtype Counts.pdf", width = 12, height = 6)
ggplot(cell_subtypes, aes(x = reorder(Cell_Type, -Cell_Count), y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "#377EB8", color = "#377EB8", width = 0.8) + 
  theme_minimal() + 
  xlab("Cell Subtype") + 
  ylab("Number of Cells") + 
  ggtitle("Cell Subtype Counts") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 20, r = 10, b = 30, l = 10, unit = "pt")
  ) +
  geom_text(aes(label = format(Cell_Count, big.mark = ",")), vjust = -0.5, color = "black", size = 5)

dev.off()

#Cell Type
cell_types <- data.frame(
  Cell_Type = c("NC", "BRAC", "CRC", "CSCC", "ESCC", "STAD", "HCC", "HNSCC", "LUAD", "LUSC", "PDAC", "PGC", "iCCA"),
  Cell_Count = c(190153, 90111, 286503, 23869, 180917, 107390, 83920, 91291, 47204, 23076, 39474, 6855, 6153)
)
cell_types$Cell_Type <- factor(cell_types$Cell_Type, 
                               levels = c("NC", "BRAC", "CRC", "CSCC", "ESCC", "STAD", "HCC", "HNSCC", "LUAD", "LUSC", "PDAC", "PGC", "iCCA"))

pdf("~/LB/data/pan cancer/CAF/fig/Fig1/Cell Type Counts.pdf", width = 12, height = 6)
ggplot(cell_types, aes(x = Cell_Type, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "#377EB8", color = "#377EB8", width = 0.8) + 
  theme_minimal() + 
  xlab("Cell Type") + 
  ylab("Number of Cells") + 
  ggtitle("Cell Type Counts") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 20, r = 10, b = 30, l = 10, unit = "pt")
  ) +
  geom_text(aes(label = format(Cell_Count, big.mark = ",")), vjust = -0.5, color = "black", size = 5) +
  scale_y_continuous(labels = scales::comma) 

dev.off()

#Tissue Type
tissue_types <- data.frame(
  Tissue_Type = c("Breast", "Colorectal", "Esophagus", "Gastric", "Larynx", "Liver", "Lung", "Lymphnode", "Oralcavity", "Oropharynx", "Peritonium", "Skin", "Pancreas"),
  Cell_Count = c(90111, 339428, 196465, 136417, 5089, 101935, 83592, 32370, 45886, 40316, 7019, 44118, 54170)
)
tissue_types$Tissue_Type <- factor(tissue_types$Tissue_Type,
                                   levels = tissue_types$Tissue_Type)

pdf("~/LB/data/pan cancer/CAF/fig/Fig1/Tissue Type Counts.pdf", width = 12, height = 6)
ggplot(tissue_types, aes(x = Tissue_Type, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "#377EB8", color = "#377EB8", width = 0.8) + 
  theme_minimal() + 
  xlab("Tissue Type") + 
  ylab("Number of Cells") + 
  ggtitle("Cell Counts by Tissue Type") + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(t = 20, r = 10, b = 30, l = 10, unit = "pt")
  ) +
  geom_text(aes(label = format(Cell_Count, big.mark = ",")), vjust = -0.5, color = "black", size = 5) +
  scale_y_continuous(labels = scales::comma) 

dev.off()




##Fig 1D##

CAF_and_Peri_harm <- qs::qread(file = '~/LB/data/pan cancer/secend data file/CAF_and_Peri_harm.qs')
CAF_harmony_samples <- readRDS(file = '~/LB/data/pan cancer/CAF/rds/CAF_harmony_samples.rds')

devtools::install_github("TheHumphreysLab/plot1cell")

Idents(CAF_and_Peri_harm)<-CAF_and_Peri_harm$megacluster
circ_data <- prepare_circlize_data(CAF_and_Peri_harm, scale = 0.6 )

cluster_colors <- rand_color(length(table(CAF_and_Peri_harm$megacluster)))
Type_colors <- rand_color(length(table(CAF_and_Peri_harm$Type)))
rep_colors <- rand_color(length(table(CAF_and_Peri_harm$orig.ident)))

colors_91 <- c("#04BF8A",'#EDAF73','#D85B5E',"#2B71A0",'#0BBDB2','#D9A2B1','#62448F','#619F94')
colors_92 <- c('#D9A2B1','#D85B5E','#62448F','#619F94',"#2B71A0","#04BF8A",'#EDAF73','#0BBDB2')
png("~/LB/data/pan cancer/CAF/fig/cluster/plot1cell_CAF_megacluster.png", width = 5, height = 5,units = 'in', res = 1000)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, 
              col.use = colors_92 ,
              bg.color = 'white', 
              kde2d.n = 1000, 
              repel = T, 
              label.cex = 1.5)

add_track(circ_data, group = "Type", colors = scPalette(length(names(table(g[["Type"]][,1])))), track_num = 2)
add_track(circ_data, group = "PatientID", colors = scPalette(length(names(table(g[["PatientID"]][,1])))), track_num = 3)
add_track(circ_data, group = "GEO", colors = scPalette(length(names(table(g[["GEO"]][,1])))), track_num = 4)
add_track(circ_data, group = "Stage", colors = scPalette(length(names(table(g[["Stage"]][,1])))), track_num = 5)

dev.off()

##Fig S2A-F##
CAF_and_Peri_harm<- RenameIdents(CAF_and_Peri_harm,`vCAF`='dCAF',`meCAF`='vCAF')
CAF_and_Peri_harm@meta.data$megacluster<-Idents(CAF_and_Peri_harm)
CAF_and_Peri_harm@meta.data$megacluster = factor(CAF_and_Peri_harm@meta.data$megacluster, levels = c('NF',"sCAF","iCAF","myCAF","vCAF","dCAF","lsCAF","Mural"))
Idents(CAF_and_Peri_harm)<-CAF_and_Peri_harm@meta.data$megacluster



pdf("~/LB/data/pan cancer/CAF/fig/cluster/umap harmony.pdf", width = 11, height = 10)
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

pdf("~/LB/data/pancancer/CAF/fig/cluster/cluster1/UMAP harmony2.pdf", width = 11, height = 10)
DimPlot(g, label = T, group.by = "rename", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = pal_50, reduction = "umap") + ggtitle("cluster")
DimPlot(g, label = T, group.by = "megacluster", label.size = 4, pt.size = 1, shuffle = T, raster = T, cols = pal_8, reduction = "umap") + ggtitle("megacluster")
dev.off()

g<-CAF_and_Peri_harm

pdf("~/LB/data/pan cancer/CAF/fig/cluster/cluster1/counts megacluster.pdf", width = 6, height = 5)
ggplot(g@meta.data, aes(x=`GEO`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`Tissue`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=Type, fill=megacluster)) + 
  geom_bar(data = filter(g@meta.data, Type != "NC"), position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis() +
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`Gender`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`Stage`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`Pathology`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`Tumor`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`T`, fill=`megacluster`)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`N`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
ggplot(g@meta.data, aes(x=`M`, fill=megacluster)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = pal_8) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
dev.off()

colors_4 <- c('#D85B5E','#EDAF73','#0BBDB2','#62448F')
pdf("~/LB/data/pan cancer/CAF/fig/cluster/cluster1/counts megacluster_stage.pdf", width = 6, height = 5)
ggplot(g@meta.data, aes(x=`megacluster`, fill=Stage)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = colors_9) + 
  theme_bw() + 
  theme(panel.grid=element_blank()) +
  RotatedAxis()+
  labs(y="Cell fraction")
dev.off()


#Fig 1F#

p1<-FeaturePlot(CAF_and_Peri_harm, features = c('PI16','CD34','HSPA6','CXCL8','CXCL14','CD74','COL11A1','MMP11','MYH11','MKI67','KRT7','RGS5'))
ggsave("~/LB/data/pan cancer/CAF/fig/cluster/CAF_marker.png",plot = p1, width = 16,height = 9,dpi = 300)



#Fig 1G#
CAF_harmony_samples <- readRDS(file = '~/LB/data/pancancer/CAF/rds/CAF_harmony_samples.rds')
CAF_harmony_samples[["RNA"]]$counts <- as(CAF_harmony_samples[["RNA"]]$counts, "dgCMatrix")
CAF_harmony_samples[["RNA"]]$data <- as(CAF_harmony_samples[["RNA"]]$data, "dgCMatrix")
write.csv(markers,file = '~/LB/data/pan cancer/CAF/marker/CAF_harmony_samples_0.3.csv')
Idents(CAF_harmony_samples) = "RNA_snn_res.0.3"
markers = FindAllMarkers(CAF_harmony_samples, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.7)
write.csv(markers,file = '~/LB/data/pan cancer/CAF/marker/CAF_harmony_samples_0.3.csv')

CAF_harmony_samples<- RenameIdents(CAF_harmony_samples,`0`='myCAF',`1`='NF',`2`='Mural',`3`='Mural',`4`='iCAF',`5`='iCAF',
                                       `7`='iCAF',`8`='NF',`9`='sCAF',`10`='iCAF',`11`='vCAF',`12`='lsCAF',`13`='iCAF',`14`='dCAF'
                                       ,`15`='lsCAF',`16`='iCAF',`17`='iCAF')

CAF_harmony_samples1 <- subset(CAF_harmony_samples,  ident = c("6"), invert = TRUE)
CAF_harmony_samples1@meta.data$megacluster<-Idents(CAF_harmony_samples1)
CAF_harmony_samples1@meta.data$megacluster = factor(CAF_harmony_samples1@meta.data$megacluster, levels = c('NF',"sCAF","iCAF","myCAF","vCAF","dCAF","lsCAF","Mural"))
saveRDS (CAF_harmony_samples1, file = '~/LB/data/pan cancer/CAF/rds/CAF_harmony_samples_rename.rds')

g<-CAF_harmony_samples1

g<- JoinLayers(g)
g <- NormalizeData(g)
g <- FindVariableFeatures(g)
g = SketchData(object = g, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
g[['sketch1']] <- 'All'
g$sketch1[colnames(g@assays$sketch)] <- 'sketch'
g1 <- subset(g, sketch1 == 'sketch')
DefaultAssay(g1) <- "RNA"
g1[['sketch']] <- NULL


######GSVA######
###### CAF全局#####
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(org.Hs.eg.db) 
library(tidyverse)
library(msigdbr)
library(GSVA)
library(limma)
library(dplyr)
library(pheatmap)
library(BiocParallel)


expr <- GetAssayData(g, slot = "data")
expr <- as.matrix(expr)

megacluster <- g@meta.data$megacluster

msgdC2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")

keggSet <- split(msgdC2$gene_symbol, msgdC2$gs_name)

zpar <- zscoreParam(expr, geneSets = keggSet)

BPPARAM <- MulticoreParam(4)

gsva_results <- gsva(zpar, verbose = TRUE)

average_scores <- rowMeans(gsva_results)

top_50_pathways <- names(sort(average_scores, decreasing = TRUE))[1:50]

top_50_scores <- gsva_results[top_50_pathways, ]

cluster_scores <- aggregate(t(top_50_scores), by = list(cluster = megacluster), FUN = mean)

heatmap_data <- as.matrix(cluster_scores[, -1])
rownames(heatmap_data) <- cluster_scores$cluster

pdf("~/LB/data/Pancancer_ICI/B_cell/fig/GSVA_KEGG_heatmap.pdf", width = 12, height = 8)
pheatmap(
  heatmap_data,
  scale = "row", 
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "GSVA Scores Heatmap"
)
dev.off()

###### iCAF亚群#####

expr <- GetAssayData(g, slot = "data")
expr <- as.matrix(expr)

rename <- g@meta.data$rename

msgdC2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")


keggSet <- split(msgdC2$gene_symbol, msgdC2$gs_name)


zpar <- zscoreParam(expr, geneSets = keggSet)


BPPARAM <- MulticoreParam(40)
BPPARAM$set(elapsedTime = 3600)
gsva_results <- gsva(zpar, BPPARAM = BPPARAM)

#gopb
msgdC5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gobpSet <- split(msgdC5$gene_symbol, msgdC5$gs_name)
gobpSet <- lapply(gobpSet, function(set) set[length(set) > 1])
gobpSet <- gobpSet[sapply(gobpSet, length) > 1]
zpar <- zscoreParam(expr, geneSets = gobpSet, minSize = 2)
zscore_results <- gsva(zpar, verbose = TRUE) 

average_scores <- rowMeans(gsva_results)

top_50_pathways <- names(sort(average_scores, decreasing = TRUE))[1:50]

top_50_scores <- gsva_results[top_50_pathways, ]

cluster_scores <- aggregate(t(top_50_scores), by = list(cluster = megacluster), FUN = mean)

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

BPPARAM <- MulticoreParam(4)

gsva_results <- gsva(zpar, BPPARAM = BPPARAM)



average_scores <- rowMeans(gsva_results)
top_50_pathways <- names(sort(average_scores, decreasing = TRUE))[1:50]

top_50_scores <- gsva_results[top_50_pathways, ]


cluster_scores <- aggregate(t(top_50_scores), by = list(cluster = rename), FUN = mean)
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


##把KEGG和GOBP放一起分析
library(GSVA)
library(msigdbr)
library(BiocParallel)
library(pheatmap)


expr <- GetAssayData(g, slot = "data")
expr <- as.matrix(expr)


megacluster <- g@meta.data$megacluster

msgdC2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
msgdBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

keggSet <- split(msgdC2$gene_symbol, msgdC2$gs_name)

gobpSet <- split(msgdBP$gene_symbol, msgdBP$gs_name)


zpar_kegg <- zscoreParam(expr, geneSets = keggSet)
zpar_gobp <- zscoreParam(expr, geneSets = gobpSet)

BPPARAM <- MulticoreParam(40)
BPPARAM$set(elapsedTime = 3600)

gsva_results_kegg <- gsva(zpar_kegg, BPPARAM = BPPARAM)
gsva_results_gobp <- gsva(zpar_gobp, BPPARAM = BPPARAM)

get_top8_per_cluster <- function(gsva_results, megacluster) {
  unique_clusters <- unique(megacluster)
  top8_per_cluster <- list()
  
  for (cluster in unique_clusters) {
    cluster_expr <- gsva_results[, megacluster == cluster]
    average_scores <- rowMeans(cluster_expr)
    top8 <- names(sort(average_scores, decreasing = TRUE))[1:8]
    top8_per_cluster[[cluster]] <- top8
  }
  
  return(top8_per_cluster)
}

top8_per_cluster_kegg <- get_top8_per_cluster(gsva_results_kegg, megacluster)
top8_per_cluster_gobp <- get_top8_per_cluster(gsva_results_gobp, megacluster)
selected_pathways_kegg <- unique(unlist(top8_per_cluster_kegg))
selected_pathways_gobp <- unique(unlist(top8_per_cluster_gobp))
selected_pathways <- c(selected_pathways_kegg, selected_pathways_gobp)

selected_scores_kegg <- gsva_results_kegg[selected_pathways_kegg, ]
selected_scores_gobp <- gsva_results_gobp[selected_pathways_gobp, ]

selected_scores <- rbind(selected_scores_kegg, selected_scores_gobp)

if (length(megacluster) != ncol(selected_scores)) {
  stop("The length of megacluster does not match the number of columns in selected_scores.")
}


cluster_scores <- aggregate(t(selected_scores), by = list(cluster = megacluster), FUN = mean)

heatmap_data <- as.matrix(cluster_scores[, -1])
rownames(heatmap_data) <- cluster_scores$cluster

pdf("~/LB/data/pan cancer/CAF/fig/GSVA/GSVA_KEGG_GOBP_heatmap_top8_per_cluster.pdf", width = 14, height = 12)
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
  main = "GSVA Scores Heatmap (Top 8 Pathways per Cluster for KEGG and GOBP)"
)
dev.off()


#GO
marker %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) -> top
ego_list = list()


for (i in unique(top$cluster)) {
  cluster_name <- as.character(i)
  
  ego_list[[cluster_name]] = enrichGO(gene = top[top$cluster == i, ]$gene, 
                                      OrgDb = "org.Hs.eg.db", 
                                      keyType = "SYMBOL", 
                                      ont = "BP")
}


ego2 <- merge_result(ego_list)
writexl::write_xlsx(ego2@compareClusterResult, "GO ORA.xlsx")
pdf("~/LB/data/pan cancer/CAF/fig/GSVA/GOBP ORA.pdf", width = 9, height = 12)
dotplot(ego2, showCategory = 10, title = "GO enrichment")
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/GSVA/GOBP ORA_8.pdf", width = 9, height = 12)
dotplot(ego2, showCategory = 8, title = "GO enrichment")
dev.off()

pdf("~/LB/data/pan cancer/CAF/fig/GSVA/GOBP ORA_5.pdf", width = 9, height = 12)
dotplot(ego2, showCategory = 5, title = "GO enrichment")
dev.off()





