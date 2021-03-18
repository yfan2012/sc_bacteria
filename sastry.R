library(tidyverse)
library(cowplot)

datadir='/dilithium/Data/NGS/projects/dunlop_rna/sastry'
dbxdir='~/Dropbox/yfan/dunlop/transcriptome'

goi=tibble(genes=c('recA', 'gadX', 'araC', 'rpoH'),
           ids=c('gene-b2699', 'gene-b3516', 'gene-b0064', 'gene-b3461'))
samps=list.dirs(file.path(datadir, 'salmon'), recursive=FALSE)

allexps=tibble(Name=as.character(),
            Length=as.numeric(),
            EffectiveLength=as.numeric(),
            TPM=as.numeric(),
            NumReads=as.numeric(),
            samp=as.character())
goiexps=tibble(Name=as.character(),
            Length=as.numeric(),
            EffectiveLength=as.numeric(),
            TPM=as.numeric(),
            NumReads=as.numeric(),
            gene=as.character(),
            samp=as.character())

for (i in samps) {
    quants=file.path(i, 'quant.sf')
    name=basename(i)

    exp=read_tsv(quants) %>%
        mutate(samp=name)
    allexps=bind_rows(allexps, exp)
    
    expgene=exp %>%
        rowwise() %>%
        filter(Name %in% goi$ids) %>%
        mutate(gene=goi$genes[goi$ids==Name])
    goiexps=bind_rows(goiexps, expgene)
}

allexps=allexps %>%
    mutate(NumReads=NumReads+1) %>%
    mutate(TPM=TPM+1)

goibarstpm=ggplot(goiexps, aes(x=gene, y=TPM, colour=samp, fill=samp, alpha=.5)) +
    geom_bar(position='dodge', stat='identity') +
    ggtitle('Genes of Interest') +
    xlab('Gene') +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
goibarsnum=ggplot(goiexps, aes(x=gene, y=NumReads, colour=samp, fill=samp, alpha=.5)) +
    geom_bar(position='dodge', stat='identity') +
    ggtitle('Genes of Interest') +
    xlab('Gene') +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
tpmdensity=ggplot(allexps, aes(x=TPM, colour=samp, fill=samp, alpha=.2)) +
    geom_density() +
    ggtitle('TPM Distributions') +
    scale_x_log10() +
    xlab('TPM') +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()
numreaddensity=ggplot(allexps, aes(x=NumReads, colour=samp, fill=samp, alpha=.2)) +
    geom_density() +
    ggtitle('Reads per gene') +
    scale_x_log10() +
    xlab('Number of Reads') +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    theme_bw()

summarypdf=file.path(dbxdir, 'sastry.pdf')    
pdf(summarypdf, h=7, w=16)
plot_grid(goibarstpm, goibarsnum)
print(tpmdensity)
print(numreaddensity)
dev.off()
