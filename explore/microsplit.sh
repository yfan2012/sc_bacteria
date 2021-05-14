#!/bin/bash

datadir=/dilithium/Data/NGS/projects/dunlop_rna/microsplit

if [ $1 == dl ] ; then
    mkdir -p $datadir/fastqs
    fastq-dump \
	-O $datadir/fastqs \
	--split-files \
	SRR11940662
fi
