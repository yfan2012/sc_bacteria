library(tidyverse)
library(tximport)
library(readr)
library(DESeq2)

datadir='/dilithium/Data/NGS/projects/dunlop_rna/deter_data'
saldir=file.path(datadir, 'salmon')
allsamps=list.dirs(saldir, recursive=FALSE)

quantfiles=file.path(allsamps, 'quant.sf')
txi=tximport(quantfiles, type='salmon', txOut=T)

