#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(MIND)
library(edgeR)

sc <- readRDS("sc.rds")
scdata<-as.matrix(sc@assays$RNA@counts)
smeta<-sc@meta.data
meta_sc<-data.frame(sample = smeta$name, cell_type = smeta$cell_type)

batch <- read.csv('stage.csv', header = TRUE)
bulk <- read.csv('tmm_matrix_combat.csv', header = TRUE)
bulk<-as.matrix(bulk)
ratio <- read.csv('ratio.csv', header = TRUE)

profile<-prior$profile
covariance<-prior$covariance
gene1<-intersect(rownames(bulk),rownames(profile))


Aest<-ratio
X<-t(bulk[match(gene1,rownames(bulk)),])
Sest<-t(profile[match(gene1,rownames(profile)),match(colnames(ratio),colnames(profile))])
lambda <- 10
iteradmm <- 10
rsCAM <- sCAMfastNonNeg(as.matrix(X), as.matrix(Aest), as.matrix(Sest), lambda = lambda, iteradmm=iteradmm, silent = T)
swcam<- rsCAM$S

rownames(swcam)<-colnames(Aest)
colnames(swcam)<-colnames(X)
dimnames(swcam)[[3]] <- rownames(X)
swcam <- aperm(swcam, perm = c(2, 1, 3))
saveRDS(swcam, file = "swcam.rds")