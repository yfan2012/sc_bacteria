#!/bin/bash

##code for looking at data from deter paper https://www.biorxiv.org/content/10.1101/2020.08.27.270272v2

datadir=/dilithium/Data/NGS/projects/dunlop_rna/deter

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

if [ $1 == aligncheck ] ; then
    for i in $datadir/align/*sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`
	outfile=$datadir/align/$prefix.stats.tsv

	echo $outfile
	
	if [ ! -s $outfile ] ; then
	    samtools flagstat \
		     -@ 24 \
		     -O tsv \
		     $i \
		     > $outfile
	fi
    done
fi

if [ $1 == alignqc ] ; then
    touch $datadir/align/alignrates.csv
    for i in $datadir/align/*stats.tsv ;
    do
	prefix=`basename $i .stats.tsv`
	rate=`head -6 $i | tail -1 | awk '{print $1}'`
	echo $prefix,$rate >> $datadir/align/alignrates.csv
    done
fi

##using reference annotation gff from ncbi and renamed to ecoli_k12.gff
if [ $1 == make_xcript_ref ] ; then
    gffread \
	-w $datadir/ref/ecoli_k12_transcriptome.fa \
	-g $ref \
	$datadir/ref/ecoli_k12.gff
	
    cp $datadir/ref/ecoli_k12_transcriptome.fa $datadir/ref/ecoli_p24KmNB82_transcriptome.fa
    grep -A 57 p24KmNB82 $datadir/ref/butzin_plasmids.fa >> $datadir/ref/ecoli_p24KmNB82_transcriptome.fa
fi

xcriptref=$datadir/ref/ecoli_p24KmNB82_transcriptome.fa

if [ $1 == salmon ] ; then
    salmon index -t $xcriptref -i $datadir/ref/ecoli_p24KmNB82_transcriptome

    mkdir -p $datadir/salmon
    
    for i in $datadir/fastqs/*_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`
	
	salmon quant \
	       -i $datadir/ref/ecoli_p24KmNB82_transcriptome \
	       -l A \
	       -1 $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
	       -2 $datadir/trimmed/${prefix}_rev_paired.fq.gz  \
	       -p 54 \
	       --validateMappings \
	       --gcBias \
	       -o $datadir/salmon/$prefix
    done
fi


if [ $1 == bowtie_xcriptome ] ; then
    bowtie2-build $xcriptref $datadir/ref/ecoli_p24KmNB82_transcriptome

    mkdir -p $datadir/align
    for i in $datadir/fastqs/SRR12518289_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`

	echo ALIGNING $prefix
	bowtie2 -p 54 \
		-x $datadir/ref/ecoli_p24KmNB82_transcriptome \
		-1 $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
		-2 $datadir/trimmed/${prefix}_rev_paired.fq.gz | \
	    samtools view -@ 54 -b | \
	    samtools sort -@ 54 -o $datadir/align/${prefix}_bowtie_xcriptome.sorted.bam
	
	samtools index $datadir/align/${prefix}_bowtie_xcriptome.sorted.bam
    done
fi


if [ $1 == align_xcriptome ] ; then
    hisat2-build $xcriptref $datadir/ref/ecoli_p24KmNB82_transcriptome

    mkdir -p $datadir/align
    for i in $datadir/fastqs/SRR12518289_1.fastq.gz ;
    do
	prefix=`basename $i _1.fastq.gz`

	echo ALIGNING $prefix
	hisat2 -p 72 \
	       -x $datadir/ref/ecoli_p24KmNB82_transcriptome \
	       -1 $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
	       -2 $datadir/trimmed/${prefix}_rev_paired.fq.gz | \
	    samtools view -@ 72 -b | \
	    samtools sort -@ 72 -o $datadir/align/${prefix}_xcriptome.sorted.bam
	
	samtools index $datadir/align/${prefix}_xcriptome.sorted.bam
    done
fi
