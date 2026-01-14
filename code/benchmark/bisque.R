#!/usr/bin/env Rscript
library(BisqueRNA)
library(Biobase)
library(Seurat)
library(data.table)
library(ggplot2)

###https://mp.weixin.qq.com/s?__biz=MzIzMjQ4ODg1Mg==&mid=2247485924&idx=1&sn=fa55d130d8457eff1e462ca94ec4961c&scene=21#wechat_redirect
###https://github.com/MedicalGenomicsLab/deconvolution_benchmarking/blob/master/src/9_2_run_bisque.R

sc<- readRDS("sc.rds") 
Idents(sc) = sc@meta.data$cell_type

sample.ids <- colnames(sc@assays$RNA@counts)
individual.labels <- sc@meta.data$orig.ident
cell.type.labels <- sc@meta.data$cell_type
sc.pheno <- data.frame(
  check.names = FALSE,
  check.rows = FALSE,
  stringsAsFactors = FALSE,
  row.names = sample.ids,
  SubjectName = individual.labels,
  cellType = cell.type.labels
)
sc.meta <- data.frame(
  labelDescription = c("SubjectName", "cellType"),
  row.names = c("SubjectName", "cellType")
)
sc.pdata <- new(
  "AnnotatedDataFrame",
  data = sc.pheno, varMetadata = sc.meta
)
scRNA.ref<-sc@assays$RNA@counts
scRNA.ref.matrix <- as.matrix(scRNA.ref)
sc.eset <- Biobase::ExpressionSet(
  assayData = scRNA.ref.matrix, phenoData = sc.pdata
)
bulk<- read.csv('bulk.csv', header = TRUE,row.names=1)
bulk.matrix <- as.matrix(bulk)
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)
gc()
res <- BisqueRNA::ReferenceBasedDecomposition(
  bulk.eset, sc.eset,
  markers = NULL, use.overlap = FALSE
)

gc()
ref.based.estimates <- res$bulk.props
write.csv(ref.based.estimates,file="ratio_bisque.csv")   
saveRDS(ref.based.estimates, file = "bulk_bisque.rds")


