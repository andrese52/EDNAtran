#!/bin/bash

#===========================================================================================================================
#GENERATING E-PROBES BASED SOLELY ON TOXIGENIC AND ATOXIGENIC STRAINS AF70 AND AF30
#===========================================================================================================================
#------------------------------------------
#1.1 data gathering
#------------------------------------------

perl genbakretriever.pl GInumbers.txt #Downloading fasta files using accession numbers
bp_genbank2gff3 --outdir . AF36.gb  ##Transforming Genbank format to GFF3, first we need to DOWNLOAD GENBANK FILE FOR THE TARGET SEQUENCE FROM NCBI
bp_genbank2gff3 --outdir . AF70.gb ##Transforming Genbank format to GFF3 AGAIN
wget ftp://ftp.tigr.org/pub/data/a_flavus/af70_20130423.assembly.fsa ##Downloading AF70 genome from JCVI
wget ftp://ftp.tigr.org/pub/data/a_flavus/af70_20130423.gff3 ##Downloading AF70 GFF3 from JCVI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000006275.2_JCVI-afl1-v2.0/GCA_000006275.2_JCVI-afl1-v2.0_genomic.fna.gz #Downloading NRRL3357 genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000952835.1_ASM95283v1/GCA_000952835.1_ASM95283v1_genomic.gff.gz #Downloading NRRL3357 gff
gunzip GCA_000006275.2_JCVI-afl1-v2.0_genomic.fna.gz #Uncompressing downloads
gunzip GCA_000952835.1_ASM95283v1_genomic.gff.gz #Uncompressing downloads
mv GCA_000006275.2_JCVI-afl1-v2.0_genomic.fna NRRL3357.fasta #Changing names
mv GCA_000952835.1_ASM95283v1_genomic.gff NRRL3357.gff #Changing names
mv af70_20130423.assembly.fsa AF70-genome.fasta #Changing names
mv af70_20130423.gff3 AF70-genome.gff3 #Changing names

#------------------------------------------
#1.2 E-probe design for the AF gene cluster
#------------------------------------------

cd /panfs/panfs.cluster/scratch/andrese52/aspergillus/input
Pipeline.pl -p AF36-20.fasta -P 20 -n AY510453.fasta -t AY510455.fasta -f AF36.gb.gff
Pipeline.pl -p AF36-40.fasta -P 40 -n AY510453.fasta -t AY510455.fasta -f AF36.gb.gff
Pipeline.pl -p AF36-60.fasta -P 60 -n AY510453.fasta -t AY510455.fasta -f AF36.gb.gff
Pipeline.pl -p AF36-80.fasta -P 80 -n AY510453.fasta -t AY510455.fasta -f AF36.gb.gff
Pipeline.pl -p AF70-20.fasta -P 20 -t AY510453.fasta -n AY510455.fasta -f AF70.gb.gff
Pipeline.pl -p AF70-40.fasta -P 40 -t AY510453.fasta -n AY510455.fasta -f AF70.gb.gff
Pipeline.pl -p AF70-60.fasta -P 60 -t AY510453.fasta -n AY510455.fasta -f AF70.gb.gff
Pipeline.pl -p AF70-80.fasta -P 80 -t AY510453.fasta -n AY510455.fasta -f AF70.gb.gff
#------------------------------------------
#1.3 E-probe design for the whole genome of AF70
#------------------------------------------


Pipeline.pl -p AF70-genome-20.fasta -P 20 -n NRRL3357.fasta -t AF70-genome.fasta -f AF70-genome.gff3 > AF70-genome-20.log
Pipeline.pl -p output/AF70-genome-40.fasta -P 40 -n input/NRRL3357.fasta -t input/AF70-genome.fasta -f input/AF70-genome.gff3
Pipeline.pl -p output/AF70-genome-60.fasta -P 60 -n input/NRRL3357.fasta -t input/AF70-genome.fasta -f input/AF70-genome.gff3
Pipeline.pl -p output/AF70-genome-80.fasta -P 80 -n input/NRRL3357.fasta -t input/AF70-genome.fasta -f input/AF70-genome.gff3
Pipeline.pl -p output/AF70-genome-100.fasta -P 100 -n input/NRRL3357.fasta -t input/AF70-genome.fasta -f input/AF70-genome.gff3


