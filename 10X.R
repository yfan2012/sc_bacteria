library(tidyverse)
library(Seurat)
library(patchwork)

dbxdir='~/Dropbox/yfan/dunlop/single_cell'
h5file=file.path(dbxdir, '10X', 'raw_feature_bc_matrix.h5')

##load h5 data
nk_all=Read10X_h5(h5file, use.names = TRUE, unique.features = TRUE)
genesums=rowSums(nk_all)
cellsums=colSums(nk_all)
starts=seq(1, dim(nk_all)[1], 6)

nk=foreach(i=starts, .combine=rbind) %dopar% {
    geneset=nk[i:i+6,]
    name=names(geneset)[1]
    gene=colSums(geneset)

##get seurat object
ecoli=CreateSeuratObject(counts=nk, project='ecoli', min.cells=3, min.features=200)



