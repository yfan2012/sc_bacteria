library(tidyverse)
library(RColorBrewer)


datadir='/dilithium/Data/NGS/projects/dunlop_rna/petriseq_data'
dbxdir='~/Dropbox/yfan/dunlop/single_cell/petriseq'

##check out supp data nk matrix, and see how many times 
exps=c('species_mix', 'growth_mix', 'growth_light_mix')
cells_per_gene=tibble(gene=as.character(),
                      numcells=as.integer(),
                      experiment=as.character(),
                      normcount=as.numeric())
for (i in exps) {
    nkfile=file.path(datadir, 'supp_data', paste0(i, '.tsv'))
    nk=read_table2(nkfile)
    nknames=c('cell', names(nk))
    nk=read_table2(nkfile, col_names=nknames, skip=1) %>%
        select(-cell)

    cells=as_tibble(enframe(colSums(nk!=0))) %>%
        mutate(value=value+1) %>%
        mutate(experiment=i) %>%
        mutate(normcount=value*1000/sum(value))
    names(cells)=c('gene', 'numcells', 'experiment', 'normcount')

    cells_per_gene=bind_rows(cells_per_gene, cells)
}

pergenefile=file.path(dbxdir, 'cells_per_gene.pdf')
pdf(pergenefile, h=7, w=12)
ggplot(cells_per_gene, aes(x=numcells, colour=experiment, fill=experiment, alpha=.25)) +
    geom_density() +
    ggtitle('Cells per gene') +
    xlab('Number of Cells') +
    scale_x_log10() +
    scale_fill_brewer(palette = 'Set2') +
    scale_colour_brewer(palette = 'Set2') +
    theme_bw()
dev.off()
