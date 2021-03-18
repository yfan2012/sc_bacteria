#!/bin/bash

datadir=/dilithium/Data/NGS/projects/dunlop_rna/sastry
ref=$datadir/ref/ecoli_k12.fa

if [ $1 == make_transcritpome ] ; then
    ##downloaded genome and annotation
    gffread \
	-w $datadir/ref/ecoli_k12_transcriptome.fa \
	-g $ref \
	$datadir/ref/ecoli_k12.gff
fi
    

if [ $1 == dl ] ; then
    ##got sra metadata table from bioproject
    ##download the k12s

    listfile=$datadir/etc/SraRunTable.txt
    
    while read p ; do
	strain=`echo $p | cut -d ',' -f 32`
	run=`echo $p | cut -d ',' -f 1`
	echo $strain
	if [ ! -z $strain ] ; then
	    if [ $strain == K-12 ] ; then
		fastq-dump \
		    -O $datadir/fastqs \
		    --split-files \
		    $run
	    fi
	fi
    done < $listfile

    pigz -p 36 $datadir/fastqs/*fastq
fi

if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed

    cp ~/software/Trimmomatic-0.39/adapters/all_adapters.fa ~/Code/yfan_ngs/dunlop_rna/
    
    for i in $datadir/fastqs/*1.fastq.gz ;
    do
	samp=`basename $i _1.fastq.gz`
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $datadir/fastqs/${samp}_1.fastq.gz $datadir/fastqs/${samp}_2.fastq.gz \
	     $datadir/trimmed/${samp}_fwd_paired.fq.gz $datadir/trimmed/${samp}_fwd_unpaired.fq.gz \
	     $datadir/trimmed/${samp}_rev_paired.fq.gz $datadir/trimmed/${samp}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done

    
fi


if [ $1 == clean ] ; then
    ##for whatever reason probs described in the paper, none of the SRR816448* reads survive trimming
    rm $datadir/trimmed/SRR816448*
fi


if [ $1 == salmon_index ] ; then
    salmon index \
	   -t $datadir/ref/ecoli_k12_transcriptome.fa \
	   -i $datadir/ref/ecoli_k12_transcriptome_salmonidx \
	   -p 36
fi


if [ $1 == salmon ] ; then
    mkdir -p $datadir/salmon
    for i in $datadir/trimmed/*_fwd_paired.fq.gz ;
    do
	samp=`basename $i _fwd_paired.fq.gz`
	salmon quant \
	      -p 36 \
	      -i $datadir/ref/ecoli_k12_transcriptome_salmonidx \
	      -l A \
	      -1 $datadir/trimmed/${samp}_fwd_paired.fq.gz \
	      -2 $datadir/trimmed/${samp}_rev_paired.fq.gz \
	      -o $datadir/salmon/$samp
    done
fi
