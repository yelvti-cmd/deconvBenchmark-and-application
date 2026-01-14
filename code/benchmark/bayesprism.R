#!/usr/bin/env Rscript
library(Seurat)
library(BayesPrism)
library(Matrix)
library(data.table)
library(ggplot2)

sc<- readRDS("sc.rds") 
Idents(sc) = sc@meta.data$cell_type
bk.dat <- read.csv('bulk.csv', header = TRUE, row.names = 1)
raw_counts <- GetAssayData(object = sc, slot = "counts")
sc.dat<-t(raw_counts)
bk.dat <-t(bk.dat)
cell.type.labels <- sc@meta.data[["cell_type"]]
cell.state.labels <- sc@meta.data[["cell_type"]]
pdf("bayesprim.pdf")
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
    sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                    species="hs", 
                                    gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","hb","act"),
                                    exp.cells=5)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
                 )
dev.off()

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,# 
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=50)
bp.res

Z.exp <- get.exp(bp = bp.res,
                        state.or.type = "type")
str(Z.exp)

saveRDS(Z.exp, file = "bulk_bayesprism.rds")

theta <- get.fraction(bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
write.csv(theta,file="ratio_bayesprism.csv")
        

