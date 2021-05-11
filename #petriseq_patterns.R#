library(tidyverse)
library(foreach)
library(doParallel)
cl=makeCluster(36)
registerDoParallel(cl, cores=36)
clusterCall(cl, function() library(tidyverse))

dbxdir='~/Dropbox/yfan/dunlop/single_cell/petriseq'
goifile=file.path(dbxdir,'growth_mix_goi.csv')

goi=read_csv(goifile)
goi[goi>1]=1
totcells=dim(goi)[1]
geneprobs=colSums(goi)/totcells


goi_collapse=goi %>%
    rowwise() %>%
    summarise(name=paste0(c_across(), collapse=''))
pop=table(goi_collapse$name)
getprob <- function(code, geneprobs) {
    ##calculates prob of each state
    ##ones where each gene has a non-zero, non-one count are a little sus
    nums=as.numeric(str_split(code, '')[[1]])
    probs=geneprobs^nums
    ones=which(probs==1)
    probs[ones]=1-geneprobs[ones]
    prob=prod(probs)
    return(prob)
}

expected=tibble(state=names(pop), numcells=pop) %>%
    rowwise() %>%
    mutate(prob=getprob(state, geneprobs)) %>%
    mutate(num_expected=prob*totcells) %>%
    mutate(ratio=round(numcells/num_expected, digits=5)) %>%
    arrange(ratio)


##sanity check with randomly generated states
x=c(0,1)
size=totcells
rando=tibble(gene1=sample(x, totcells, replace=TRUE, prob=c(1-geneprobs[1], geneprobs[1]))) %>%
    mutate(gene2=sample(x, totcells, replace=TRUE, prob=c(1-geneprobs[2], geneprobs[2]))) %>%
    mutate(gene3=sample(x, totcells, replace=TRUE, prob=c(1-geneprobs[3], geneprobs[3]))) %>%
    mutate(gene4=sample(x, totcells, replace=TRUE, prob=c(1-geneprobs[4], geneprobs[4]))) %>%
    mutate(gene5=sample(x, totcells, replace=TRUE, prob=c(1-geneprobs[5], geneprobs[5]))) %>%
    mutate(gene6=sample(x, totcells, replace=TRUE, prob=c(1-geneprobs[6], geneprobs[6])))
randoprobs=colSums(rando)/totcells
rando=rando %>%
    rowwise() %>%
    summarise(name=paste0(c_across(), collapse=''))
randopop=table(rando$name)
randoexpected=tibble(state=names(randopop), numcells=randopop) %>%
    rowwise() %>%
    mutate(prob=getprob(state, randoprobs)) %>%
    mutate(num_expected=prob*totcells) %>%
    mutate(ratio=round(numcells/num_expected, digits=5)) %>%
    arrange(ratio)


##try permutation test?
permute_genes <- function(goi) {
    goipermute=goi %>%
        mutate(across(everything(), sample)) %>%
        rowwise() %>%
        summarise(name=paste0(c_across(), collapse=''))
    permutepop=table(goipermute$name)
    permute_expected=tibble(state=names(permutepop), numcells=permutepop) %>%
        rowwise() %>%
        mutate(prob=getprob(state, geneprobs)) %>%
        mutate(num_expected=prob*totcells) %>%
        mutate(ratio=round(numcells/num_expected, digits=5)) %>%
        arrange(ratio)
    return(permute_expected)
}

library(binaryLogic)

info=tibble(nums=seq(0, 63, 1)) %>%
    rowwise() %>%
    mutate(states=paste(as.character(as.binary(nums, n=6)), collapse=''))

perms=foreach(i=1:50000, .combine=cbind) %dopar% {
    perm=permute_genes(goi)
    result=matrix(0, nrow=64, ncol=1)
    rownames(result)=info$states
    result[perm$state,1]=perm$ratio
    return(result)
}
    
toplot=tibble('g0'=perms['000000',],
              'g1'=perms['000001',],
              'g2'=perms['000010',],
              'g3'=perms['000100',],
              'g4'=perms['001000',],
              'g5'=perms['010000',],
              'g6'=perms['100000',],
              'c7'=perms['001100',])



library(RColorBrewer)
library(cowplot)
plot_permute <- function(toplotgene, state) {
    geneplot=tibble(gene=toplotgene)

    bound=expected$ratio[expected$state==state]
    pval=sum(geneplot>bound)/dim(geneplot)[1]
    
    plot=ggplot(geneplot, aes(x=gene, alpha=.3, colour='#66C2A5', fill='#66C2A5')) +
        geom_histogram() +
        ggtitle(paste0('expression state: ', state)) +
        geom_vline(xintercept=bound) +
        xlim(.7,1.2) +
        scale_colour_brewer(palette='Set2') +
        scale_fill_brewer(palette='Set2') +
        theme_bw()
    return(plot)
}

g1plot=plot_permute(toplot$g1, '000001')
g2plot=plot_permute(toplot$g2, '000010')
g3plot=plot_permute(toplot$g3, '000100')
g4plot=plot_permute(toplot$g4, '001000')
g5plot=plot_permute(toplot$g5, '010000')
g6plot=plot_permute(toplot$g6, '100000')



permpdf=file.path(dbxdir,'perm_distributions.pdf')
pdf(permpdf, h=6, w=17)
plot_grid(g1plot, g2plot, g3plot, g4plot, g5plot, g6plot, ncol=3, align='v')
dev.off()




        



    
    






