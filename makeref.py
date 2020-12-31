from Bio import SeqIO
import pysam

refdir='/dilithium/Data/NGS/projects/dunlop_rna/deter_data/ref'
ecolifile=refdir+'/ecoli_k12.fa'
plasgbfile=refdir+'/Zip_1_Plasmids.gb'
plasfile=refdir+'/butzin_plasmids.fa'

##had to manually add tabs before length field in each locus line to get the right format
SeqIO.convert(plasgbfile, "genbank", plasfile, "fasta")

plas=pysam.FastaFile(plasfile)
seqs=plas.references
plasdict={ x:plas.fetch(x) for x in seqs }
ecoli=pysam.FastaFile(ecolifile)
chrs=ecoli.references
ecolidict={ x:ecoli.fetch(x) for x in chrs }

deterref=refdir+'/ecoli_p24KmNB82.fa'
with open(deterref, 'w') as f:
    for i in ecolidict:
        f.write('>'+i+'\n')
        f.write(ecolidict[i]+'\n')
    f.write('>p24KmNB82'+'\n')
    f.write(plasdict['p24KmNB82']+'\n')
    f.close()
    
