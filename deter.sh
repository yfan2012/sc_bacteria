#!/bin/bash

##code for looking at data from deter paper https://www.biorxiv.org/content/10.1101/2020.08.27.270272v2

datadir=/dilithium/Data/NGS/projects/dunlop_rna/deter_data

if [ $1 == dl ] ; then
    listfile=$datadir/accession_list.txt

    while read p ; do
	fastq-dump \
	    -O $datadir/fastqs \
	    --split-files \
	    $p
    done < $listfile

fi

if [ $1 == zip ] ; then
    pigz -p 54 $datadir/fastqs/*fastq
fi

if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed
    cp ~/software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa ~/Code/yfan_ngs/dunlop_rna/
    
    for i in $datadir/fastqs/*_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`
	
        java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
             -threads 36 -phred33 \
             $datadir/fastqs/${prefix}_1.fastq.gz $datadir/fastqs/${prefix}_2.fastq.gz \
             $datadir/trimmed/${prefix}_fwd_paired.fq.gz $datadir/trimmed/${prefix}_fwd_unpaired.fq.gz \
             $datadir/trimmed/${prefix}_rev_paired.fq.gz $datadir/trimmed/${prefix}_rev_unpaired.fq.gz \
             ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == makeref ] ; then
    python makeref.py
fi

ref=$datadir/ref/ecoli_p24KmNB82.fa
deterref=$datadir/ref/CFP-LAA-genome.fa

if [ $1 == align ] ; then
    hisat2-build $ref $datadir/ref/ecoli_p24KmNB82

    mkdir -p $datadir/align
    for i in $datadir/fastqs/*_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`

	echo ALIGNING $prefix
	hisat2 -p 72 \
	       -x $datadir/ref/ecoli_p24KmNB82 \
	       -1 $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
	       -2 $datadir/trimmed/${prefix}_rev_paired.fq.gz | \
	    samtools view -@ 72 -b | \
	    samtools sort -@ 72 -o $datadir/align/$prefix.sorted.bam
	
	samtools index $datadir/align/$prefix.sorted.bam
    done
fi

if [ $1 == bowtie ] ; then
    bowtie2-build $ref $datadir/ref/ecoli_p24KmNB82

    mkdir -p $datadir/align
    for i in $datadir/fastqs/SRR12518289_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`

	echo ALIGNING $prefix
	bowtie2 -p 54 \
		-x $datadir/ref/ecoli_p24KmNB82 \
		-1 $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
		-2 $datadir/trimmed/${prefix}_rev_paired.fq.gz | \
	    samtools view -@ 54 -b | \
	    samtools sort -@ 54 -o $datadir/align/${prefix}_bowtie.sorted.bam
	
	samtools index $datadir/align/${prefix}_bowtie.sorted.bam

    done
fi

if [ $1 == align_deterref ] ; then
    hisat2-build $deterref $datadir/ref/CFP-LAA-genome

    mkdir -p $datadir/align
    for i in $datadir/fastqs/SRR12518289_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`

	echo ALIGNING $prefix
	hisat2 -p 72 \
	       -x $datadir/ref/ecoli_p24KmNB82 \
	       -1 $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
	       -2 $datadir/trimmed/${prefix}_rev_paired.fq.gz | \
	    samtools view -@ 72 -b | \
	    samtools sort -@ 72 -o $datadir/align/${prefix}_deterref.sorted.bam
	
	samtools index $datadir/align/${prefix}_deterref.sorted.bam
    done
fi

if [ $1 == bowtie_deterref ] ; then
    bowtie2-build $deterref $datadir/ref/CFP-LAA-genome

    mkdir -p $datadir/align
    for i in $datadir/fastqs/SRR12518289_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`

	echo ALIGNING $prefix
	bowtie2 -p 54 \
		-x $datadir/ref/CFP-LAA-genome \
		-1 $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
		-2 $datadir/trimmed/${prefix}_rev_paired.fq.gz | \
	    samtools view -@ 54 -b | \
	    samtools sort -@ 54 -o $datadir/align/${prefix}_bowtie_deterref.sorted.bam
	
	samtools index $datadir/align/${prefix}_bowtie_deterref.sorted.bam

    done
fi

