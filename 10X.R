library(tidyverse)
library(Seurat)

dbxdir='~/Dropbox/yfan/dunlop/single_cell'
h5file=file.path(dbxdir, '10X', 'raw_feature_bc_matrix.h5')

test=Read10X_h5(h5file, use.names = TRUE, unique.features = TRUE)

