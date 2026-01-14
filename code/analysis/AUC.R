#!/usr/bin/env Rscript
#conda install -c bioconda bioconductor-aucell=1.20.1 
library(AUCell)
library(Seurat)
library(ggplot2)
library(dplyr)

sc <- readRDS("sc.rds")
Idents(sc)=sc@meta.data$cellType

# 1. 
# module
module <- read.csv('wgcna_module.csv', header = TRUE)
geneSets <- split(module$X, module$module)

# 2.
exprMatrix <- GetAssayData(sc, assay = "RNA", slot = "data")  # log-normalized
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)

# 3. AUC
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# 4. Seurat
auc_matrix <- t(as.matrix(getAUC(cells_AUC)))
sc@meta.data <- cbind(sc@meta.data, auc_matrix)

# 5. 
modules <- unique(module$module)

pdf("UMAP_AUC.pdf")
plot_list <- list()
for (mod in modules) {
  p <- FeaturePlot(sc, features = mod, reduction = "umap", raster = TRUE) +
    scale_colour_gradientn(colours = c("blue", "white", "red")) +
    ggtitle(paste("Module:", mod)) +
    theme(plot.title = element_text(hjust = 0.5))
  plot_list[[mod]] <- p
}
CombinePlots(plot_list)
dev.off()

# 6. 
auc_matrix_df <- data.frame(auc_matrix)
auc_matrix_df$cellType <- sc@meta.data$cellType
write.csv(auc_matrix_df, file = "module_activity_scores.csv")