#!/usr/bin/env Rscript

library(Seurat)
library(MuSiC)
library(Biobase)
library(SingleCellExperiment)
library(scater)

##https://xuranw.github.io/MuSiC/articles/MuSiC.html
###https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7181686/
sc<- readRDS("sc.rds") 
Idents(sc) = sc@meta.data$cell_type
table(sc@meta.data$cell_type)
sc <- as.SingleCellExperiment(sc)

bk.dat <- read.csv('bulk.csv', header = TRUE,row.names=1)
object <- new("ExpressionSet", exprs=as.matrix(bk.dat))
bk.dat <- exprs(object)

prop = music_prop(bulk.mtx = bk.dat, sc.sce = sc, clusters = 'cell_type',
                               samples = 'name', verbose = F)$Est.prop.weighted
prop_all = cbind('proportion'=c(prop),
                 'celltype'=rep(colnames(prop), 
                                each=nrow(prop)), 
                 'sampleID'=rep(rownames(prop),times=ncol(prop)),
                 'Method'='MuSiC')
prop_all=as.data.frame(prop_all)
prop_all$proportion=as.numeric(as.character(prop_all$proportion))
write.csv(prop_all,file="ratio_MuSiC.csv")
