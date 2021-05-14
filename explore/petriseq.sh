#!/bin/bash

datadir=/dilithium/Data/NGS/projects/dunlop_rna/petriseq
dbxdir=~/Dropbox/yfan/dunlop/single_cell/petriseq

if [ $1 == grab_supp ] ; then
    wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-020-0729-6/MediaObjects/41564_2020_729_MOESM4_ESM.gz
fi

if [ $1 == sep_experiments ] ; then
    ##separate out datasets according to supp table 6 information
    head -1 $datadir/supp_data/41564_2020_729_MOESM4_ESM > $datadir/supp_data/species_mix.tsv
    grep SB346 $datadir/supp_data/41564_2020_729_MOESM4_ESM >> $datadir/supp_data/species_mix.tsv
    
    head -1 $datadir/supp_data/41564_2020_729_MOESM4_ESM > $datadir/supp_data/growth_light_mix.tsv
    grep 394A $datadir/supp_data/41564_2020_729_MOESM4_ESM >> $datadir/supp_data/growth_light_mix.tsv

    head -1 $datadir/supp_data/41564_2020_729_MOESM4_ESM > $datadir/supp_data/growth_mix.tsv
    grep SB442 $datadir/supp_data/41564_2020_729_MOESM4_ESM >> $datadir/supp_data/growth_mix.tsv
fi

if [ $1 == dl_sra ] ; then
    listfile=$dbxdir/sra_accessionlist.txt

    while read p ; do
	fastq-dump \
	    -O $datadir/raw \
	    --split-files \
	    $p
    done < $listfile
fi


