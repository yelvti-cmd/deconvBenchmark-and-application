#!/usr/bin/env Rscript

library(Seurat)
library(MuSiC)
library(Biobase)
library(SingleCellExperiment)
library(scater)

##https://xuranw.github.io/MuSiC/articles/MuSiC.html
###https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7181686/
sc<- readRDS("/public/home/2018013/02users/yewen/mammary_harmony.rds") 
cell_type <- c("T cells", "LumSec", "Basal cells", "T cells", "Myoepithelial cells","Blood vascular endothelial cells", "Macrophages","LumHR","Basal cells", "Dendritic cells","Fibroblasts", "T cells",
                 "Neutrophils", "Lymphatic endothelial cells", "Myoepithelial cells", "Proliferative T cells", "LumSec", "Mast cells","Basal cells", "T cells", "LumSec", "LumSec")

names(cell_type) <- levels(sc)
sc <- RenameIdents(sc, cell_type)
sc@meta.data$cell_type <- Idents(sc)
Idents(sc) = sc@meta.data$cell_type
table(sc@meta.data$cell_type)

sc <- as.SingleCellExperiment(sc)

bk.dat <- read.csv('/public/home/2018013/02users/yewen/merged_gene_count_matrix_bulk_add.csv', header = TRUE,row.names=1)
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
write.csv(prop_all,file="/public/home/2018013/02users/yewen/ratio_MuSiC.csv")

pdf("ratio_MuSiC.pdf")
ggplot(prop_all) + 
  geom_bar(aes(x = sampleID,y = proportion,fill = celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()