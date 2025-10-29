

#####所有的数据_72samples#####

BC_all <- read.csv("~/LB/data/CODEX/data/BC_all.csv", header = TRUE)
CRC_all <- read.csv("~/LB/data/CODEX/data/CRC_all.csv", header = TRUE)
CSCC_all <- read.csv("~/LB/data/CODEX/data/CSCC_all.csv", header = TRUE)
GC_all_1 <- read.csv("~/LB/data/CODEX/data/GC_all_1.csv", header = TRUE)
GC_all_2 <- read.csv("~/LB/data/CODEX/data/GC_all_2.csv", header = TRUE)
HCC_all <- read.csv("~/LB/data/CODEX/data/HCC_all.csv", header = TRUE)
LC_all <- read.csv("~/LB/data/CODEX/data/LC_all.csv", header = TRUE)
EC_all_1 <- read.csv("~/LB/data/CODEX/data/EC_all_1.csv", header = TRUE)
HNSCC_all <- read.csv("~/LB/data/CODEX/data/HNSCC_all.csv", header = TRUE)

merged_data <- rbind(BC_all, CRC_all, CSCC_all, GC_all_1, GC_all_2, HCC_all, LC_all,EC_all_1,HNSCC_all )



merged_data <- merged_data %>%
  mutate(Classification = str_replace_all(Classification, "\\s*\\(.*?\\)", ""))


mapping <- data.frame(
  Classification = c(
    "Pan-Cytokeratin",
    "CD31",
    "CD68",
    "iNOS: CD68",
    "CD163",
    "CD14: MPO",
    "CD14",
    "CD4",
    "CD8",
    "FOXP3: CD4",
    "CD8: Granzyme B",
    "CD8: CD45RO",
    "CD45RO: CD4","CD8: TOX",
    "CD79a",
    "HIF1A",
    "IFNG",
    "MMP11: SMA",
    "SMA: Collagen IV",
    "Ki67: SMA",
    "Keratin 5: Pan-Cytokeratin: SMA",
    "GM-CSF: MCP-2: SMA",
    "SMA: CXCL14",
    "GM-CSF: SMA",
    "G-CSF: SMA",
    "SMA: IDO1",
    "SMA: HLA-DR"
  ),
  Name = c(
    "Epithelial cells",
    "Endothelial cells",
    "macrophages",
    "M1 macrophages",
    "M2 macrophages",
    "Neutrophils",
    "Monocytes",
    "CD4 T cells",
    "CD8 T cells",
    "Treg",
    "Cytotoxic T cells",
    "Memory T cells","Memory T cells",
    "Tex",
    "B cells",
    "Hypoxia",
    "IFNG",
    "MyCAF",
    "vCAF",
    "dCAF",
    "lsCAF",
    "Mono_iCAF",
    "Mac1_iCAF",
    "Mac2_iCAF",
    "Gran_iCAF",
    "Tem_iCAF",
    "apCAF"
  )
)

# 初始化一个空的数据框来存储结果
result <- data.frame()

# 遍历映射关系表，筛选并添加名称
for (i in 1:nrow(mapping)) {
  temp <- merged_data %>%
    filter(str_detect(`Classification`, sprintf("^%s$", mapping$Classification[i]))) %>%
    mutate(Name = mapping$Name[i])
  result <- bind_rows(result, temp)
}

write.csv(result, '~/LB/data/CODEX/data/all_72_select.csv', row.names = FALSE)


#######第一张片子_36samples#####
merged_data <- rbind(BC_all, GC_all_1,LC_all,EC_all_1 )



merged_data <- merged_data %>%
  mutate(Classification = str_replace_all(Classification, "\\s*\\(.*?\\)", ""))


mapping <- data.frame(
  Classification = c(
    "Pan-Cytokeratin",
    "CD31",
    "CD68",
    "iNOS: CD68",
    "CD163",
    "CD14: MPO",
    "CD14",
    "CD4",
    "CD8",
    "FOXP3: CD4",
    "CD8: Granzyme B",
    "CD8: CD45RO",
    "CD45RO: CD4","CD8: TOX",
    "CD79a",
    "HIF1A",
    "IFNG",
    "MMP11: SMA",
    "SMA: Collagen IV",
    "Ki67: SMA",
    "Keratin 5: Pan-Cytokeratin: SMA",
    "GM-CSF: MCP-2: SMA",
    "SMA: CXCL14",
    "GM-CSF: SMA",
    "G-CSF: SMA",
    "SMA: IDO1",
    "SMA: HLA-DR"
  ),
  Name = c(
    "Epithelial cells",
    "Endothelial cells",
    "macrophages",
    "M1 macrophages",
    "M2 macrophages",
    "Neutrophils",
    "Monocytes",
    "CD4 T cells",
    "CD8 T cells",
    "Treg",
    "Cytotoxic T cells",
    "Memory T cells","Memory T cells",
    "Tex",
    "B cells",
    "Hypoxia",
    "IFNG",
    "MyCAF",
    "vCAF",
    "dCAF",
    "lsCAF",
    "Mono_iCAF",
    "Mac1_iCAF",
    "Mac2_iCAF",
    "Gran_iCAF",
    "Tem_iCAF",
    "apCAF"
  )
)

# 初始化一个空的数据框来存储结果
result <- data.frame()

# 遍历映射关系表，筛选并添加名称
for (i in 1:nrow(mapping)) {
  temp <- merged_data %>%
    filter(str_detect(`Classification`, sprintf("^%s$", mapping$Classification[i]))) %>%
    mutate(Name = mapping$Name[i])
  result <- bind_rows(result, temp)
}

write.csv(result, '~/LB/data/CODEX/data/all_36_1_select.csv', row.names = FALSE)



########第二张片子_36sampeles#######
merged_data <- rbind(HNSCC_all, CRC_all, CSCC_all, GC_all_2, HCC_all)



merged_data <- merged_data %>%
  mutate(Classification = str_replace_all(Classification, "\\s*\\(.*?\\)", ""))


mapping <- data.frame(
  Classification = c(
    "Pan-Cytokeratin",
    "CD31",
    "CD68",
    "iNOS: CD68",
    "CD163",
    "CD14: MPO",
    "CD14",
    "CD4",
    "CD8",
    "FOXP3: CD4",
    "CD8: Granzyme B",
    "CD8: CD45RO",
    "CD45RO: CD4","CD8: TOX",
    "CD79a",
    "HIF1A",
    "IFNG",
    "MMP11: SMA",
    "SMA: Collagen IV",
    "Ki67: SMA",
    "Keratin 5: Pan-Cytokeratin: SMA",
    "GM-CSF: MCP-2: SMA",
    "SMA: CXCL14",
    "GM-CSF: SMA",
    "G-CSF: SMA",
    "SMA: IDO1",
    "SMA: HLA-DR"
  ),
  Name = c(
    "Epithelial cells",
    "Endothelial cells",
    "macrophages",
    "M1 macrophages",
    "M2 macrophages",
    "Neutrophils",
    "Monocytes",
    "CD4 T cells",
    "CD8 T cells",
    "Treg",
    "Cytotoxic T cells",
    "Memory T cells","Memory T cells",
    "Tex",
    "B cells",
    "Hypoxia",
    "IFNG",
    "MyCAF",
    "vCAF",
    "dCAF",
    "lsCAF",
    "Mono_iCAF",
    "Mac1_iCAF",
    "Mac2_iCAF",
    "Gran_iCAF",
    "Tem_iCAF",
    "apCAF"
  )
)

# 初始化一个空的数据框来存储结果
result <- data.frame()

# 遍历映射关系表，筛选并添加名称
for (i in 1:nrow(mapping)) {
  temp <- merged_data %>%
    filter(str_detect(`Classification`, sprintf("^%s$", mapping$Classification[i]))) %>%
    mutate(Name = mapping$Name[i])
  result <- bind_rows(result, temp)
}

write.csv(result, '~/LB/data/CODEX/data/all_36_2_select.csv', row.names = FALSE)





#Fig 3C

data <- read.csv("~/LB/data/CODEX/data/heatmap/heatmap_72_all_marker_select.csv", row.names = 1, header = TRUE)

# 将数据转换为数值矩阵（如果需要）
data_matrix <- as.matrix(data)


pdf("~/LB/data/CODEX/fig/all/all_72_heatmap_cell_marker.pdf", width = 10, height = 8)
pheatmap(data_matrix,
         scale = "none",  # 不进行标准化
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 10,
         border_color = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 15,    # 设定每个热块的宽度
         cellheight = 15,   # 设定每个热块的高度
         col = colorRampPalette(c( "white", "#F29F7C","#860422"))(100),
         breaks =seq(0, 0.6, length.out = 100)
)
dev.off()



#Fig 3D

data <- read.csv("~/LB/data/CODEX/data/heatmap/EC_all1.csv", row.names = 1, header = TRUE)

# 将数据转换为数值矩阵（如果需要）
data_matrix <- as.matrix(data)

pdf("~/LB/data/CODEX/fig/all/all_72_heatmap_cell.pdf", width = 10, height = 8)
pheatmap(data_matrix,
         scale = "none",  # 不进行标准化
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 10,
         border_color = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 15,    # 设定每个热块的宽度
         cellheight = 15,   # 设定每个热块的高度
         col = colorRampPalette(c( "white", "#F29F7C","#860422"))(100),
         breaks =seq(0, 1, length.out = 100)
)
dev.off()



























































