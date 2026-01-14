#!/usr/bin/env Rscript
#conda install -c bioconda bioconductor-aucell=1.20.1  # 对应 R 4.2.x
# 加载必要的包
library(AUCell)
library(Seurat)
library(ggplot2)
library(dplyr)

sc <- readRDS("/public/home/2018013/02users/yewen/CATD_snakemake/Input/Cell_splits/sc_C0_seeded.rds")
Idents(sc)=sc@meta.data$cellType
sc<-subset(sc, subset = cellType != "Prolifer" )
sc@meta.data$cellType[which(sc@meta.data$cellType == "NKT")] <- "T"
sc@meta.data$cellType <- droplevels(sc@meta.data$cellType)
# 1. 准备基因集
# 从module数据框中提取基因集
module <- read.csv('/public/home/2018013/02users/yewen/beifen/wgcna_keymodule_422.csv', header = TRUE)
geneSets <- split(module$X, module$module)

# 2. 计算基因排名
# 使用表达矩阵（假设已经normalized和log transformed）
exprMatrix <- GetAssayData(sc, assay = "RNA", slot = "data")  # 使用log-normalized数据
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)

# 3. 计算AUC值
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# 4. 将结果添加到Seurat对象中
auc_matrix <- t(as.matrix(getAUC(cells_AUC)))
sc@meta.data <- cbind(sc@meta.data, auc_matrix)

# 5. 可视化每个模块的活性
# 假设你已经运行过UMAP降维（如果没有需要先运行RunUMAP）

# 获取模块名称
modules <- unique(module$module)

# 为每个模块绘制UMAP图
pdf("UMAP_AUC.pdf")
plot_list <- list()
for (mod in modules) {
  p <- FeaturePlot(sc, features = mod, reduction = "umap", raster = TRUE) +
    scale_colour_gradientn(colours = c("blue", "white", "red")) +
    ggtitle(paste("Module:", mod)) +
    theme(plot.title = element_text(hjust = 0.5))
  plot_list[[mod]] <- p
}

# 显示所有模块的UMAP图
CombinePlots(plot_list)
dev.off()
# 6. (可选) 保存AUC矩阵
auc_matrix_df <- data.frame(auc_matrix)
auc_matrix_df$cellType <- sc@meta.data$cellType
write.csv(auc_matrix_df, file = "module_activity_scores.csv")