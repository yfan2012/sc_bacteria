#!/bin/bash

##code for looking at data from deter paper https://www.biorxiv.org/content/10.1101/2020.08.27.270272v2

datadir=/dilithium/Data/NGS/projects/dunlop_rna/deter_data

if [ $1 == dl ] ; then
    listfile=$datadir/accession_list.txt

    mkdir -p $datadir/fastqs
    while read p ; do
	fastq-dump --outdir $datadir/fastqs $p
    done < $listfile
fi

if [ $1 == align ] ; then
    hisat2-build $datadir/ref/ecoli_k12.fa $datadir/ref/CFP-LAA-genome

    mkdir -p $datadir/align
    for i in $datadir/fastqs/*fastq ;
    do
	prefix=`basename $i .fastq`

	if [ ! -f $prefix.sorted.bam.bai ] ; then

	    echo ALIGNING $prefix
	    hisat2 -p 36 \
		   -x $datadir/ref/CFP-LAA-genome \
		   -U $i | \
		samtools view -@ 36 -b | \
		samtools sort -@ 36 -o $datadir/align/$prefix.sorted.bam

	    samtools index $datadir/align/$prefix.sorted.bam
	fi
    done
fi

if [ $1 == bowtie ] ; then
    bowtie2-build $datadir/ref/ecoli_k12.fa $datadir/ref/CFP-LAA-genome

    mkdir -p $datadir/align
    for i in $datadir/fastqs/*fastq ;
    do
	prefix=`basename $i .fastq`

	if [ ! -f $prefix.sorted.bam.bai ] ; then

	    echo ALIGNING $prefix
	    bowtie2 -p 24 \
		   -x $datadir/ref/CFP-LAA-genome \
		   -U $i | \
		samtools view -@ 24 -b | \
		samtools sort -@ 24 -o $datadir/align/$prefix.sorted.bam

	    samtools index $datadir/align/$prefix.sorted.bam
	fi
    done
fi

