#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(MIND)
library(edgeR)

sc <- readRDS("sc.rds")
scdata<-as.matrix(sc@assays$RNA@counts)
smeta<-sc@meta.data
meta_sc<-data.frame(sample = smeta$name, cell_type = smeta$cell_type)

prior = get_prior(sc = scdata, meta_sc = meta_sc,filter_pd = T)

batch <- read.csv('stage.csv', header = TRUE)
bulk <- read.csv('bulk.csv', header = TRUE)
bulk<-as.matrix(bulk)
ratio <- read.csv('ratio.csv', header = TRUE)  ###sample x cell type

profile<-prior$profile
covariance<-prior$covariance
gene1<-intersect(rownames(bulk),rownames(profile))
bmind = bMIND(bulk = bulk[match(gene1,rownames(bulk)),], frac = as.matrix(ratio), sample_id = rownames(ratio), y=batch[,2], profile=profile[match(gene1,rownames(profile)),match(colnames(ratio),colnames(profile))],  covariance = covariance[match(gene1,rownames(covariance)),match(colnames(ratio),colnames(profile)),match(colnames(ratio),colnames(profile))],ncore = 30)
saveRDS(bmind, file = "bmind.rds")