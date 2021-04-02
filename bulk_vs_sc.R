library(tidyverse)
library(cowplot)

datadir='/dilithium/Data/NGS/projects/dunlop_rna'
dbxdir='~/Dropbox/yfan/dunlop'

##read sc info
scinfo=read_csv(file.path(dbxdir, 'single_cell', 'petriseq', 'sc_rank.csv')) %>%
    rowwise() %>%
    mutate(gene=str_split(name, ':')[[1]][2])

##read bulk info
bulkinfo=read_csv(file.path(dbxdir, 'transcriptome', 'goirank.csv')) %>%
    rename(experiment=samp) %>%
    mutate(tech='bulk') %>%
    select(gene, experiment, rank, tech)

genekey=tibble(short=unique(bulkinfo$gene)[order(unique(bulkinfo$gene))],
               long=unique(scinfo$gene)[order(unique(scinfo$gene))])

scinfo=scinfo %>%
    rowwise() %>%
    mutate(gene=genekey$short[genekey$long==gene]) %>%
    mutate(tech='sc') %>%
    select(gene, experiment, rank, tech)

both=bind_rows(scinfo, bulkinfo)


rankpdf=file.path(dbxdir, 'goi_bulk_vs_sc.pdf')
pdf(rankpdf, h=7, w=12)
ggplot(both, aes(x=gene, y=rank, colour=experiment, fill=experiment, alpha=.3)) +
    geom_bar(position='dodge', stat='identity') +
    ggtitle('GOI ranks') +
    ##scale_y_reverse() +
    scale_fill_brewer(palette='Set2') +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
dev.off()


scexps=c('growth_light_mix', 'species_mix', 'growth_mix')
for (i in scexps) {
    csv=file.path(dbxdir, 'single_cell', 'petriseq', paste0(i, '.genes.csv'))
    genes=read_csv(csv) %>%
        arrange(-avg)
    write_csv(genes, csv)
}
