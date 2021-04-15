library(tidyverse)

goifile='~/Dropbox/yfan/dunlop/single_cell/petriseq/growth_mix_goi.csv'

geneprobs=colSums(goi)/numcells
numcells=dim(goi)[1]

getprob <- function(code, geneprobs) {
    nums=as.numeric(str_split(code, '')[[1]])
    probs=geneprobs^nums
    ones=which(probs==1)
    probs[ones]=1-geneprobs[ones]
    prob=prod(probs)
    return(prob)
}

goi=read_csv(goifile) %>%
    mutate_all(as.character) %>%
    rowwise() %>%
    summarise(name=paste0(c_across(), collapse=''))
pop=table(goi$name)
expected=tibble(names=names(pop), pop=pop) %>%
    rowwise() %>%
    mutate(prob=getprob(names, geneprobs)) %>%
    mutate(num=prob*numcells) %>%
    mutate(ratio=log(pop/num, base=10))
    
    





        



    
    






