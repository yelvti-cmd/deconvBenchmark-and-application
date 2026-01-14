#!/usr/bin/env Rscript
library(Seurat)
library(preprocessCore)
library(scAB)
library(Scissor)
library(gplots)
#options(future.globals.maxSize = 400 * 1024^3)

#### =============================================
#### 1. read data
#### =============================================
#### 1.1 scRNA
sc_dataset <- readRDS("sc.rds")
Idents(sc_dataset) <- 'cell_type'
meta_data <- sc@meta.data

sc_dataset@graphs <- list()
sc_dataset@commands <- list()
UMAP_celltype <- DimPlot(sc_dataset, reduction ="umap", group.by="cell_type",label = T)
sc_dataset@graphs <- list()
sc_dataset@commands <- list()
sce <- run_seurat(sc_dataset,verbose = FALSE)

#### 1.2 bulk RNA
bulk <- read.csv('bulk.csv', header = TRUE)
bulk_dataset <- as.matrix(bulk)
#### 1.3 phenotype
phenotype <- read.csv('stage.csv', header = TRUE)
phenotype <- phenotype[,"sort"]
#dry=0, lactating=1
tag <- c("dry","lactating")
gc()

#### =============================================
#### 2. Scissor
#### =============================================
source('scissor_function.R')  
infos1 <- Scissor(bulk_dataset, sce, phenotype, 
                     tag = tag,
                     alpha = 0.05, 
                     cutoff = 0.03, #the number of the Scissor selected cells should not exceed 20% of total cells in the single-cell data
                     family = "binomial", 
                     Save_file = './Scissor.RData',
                     Load_file = NULL
)
gc()
Scissor_select <- rep(0, ncol(sce))
names(Scissor_select) <- colnames(sce)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"
sce <- AddMetaData(sce, metadata = Scissor_select, col.name = "scissor")

UMAP_seurat <- DimPlot(sce, reduction ="umap", group.by="cell_type",label = T)
UMAP_scissor <- DimPlot(sce, reduction = 'umap', group.by = 'scissor', cols = c('grey','royalblue','indianred1'), pt.size = 0.001, order = c("Scissor+","Scissor-"))
cairo_pdf('Scissor_umap.pdf', width=15)
patchwork::wrap_plots(plots = list(UMAP_celltype, UMAP_seurat, UMAP_scissor), ncol = 3)
dev.off()
cairo_pdf('Scissor_stat.pdf', width=15)
balloonplot(table(sce$scissor,sce$cell_type))
dev.off()