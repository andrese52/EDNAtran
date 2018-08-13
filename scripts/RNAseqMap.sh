#!/bin/bash

#-------------------------------------------------------------------------------------------------------
#Most of these is taken from here: http://www.bioconductor.org/help/workflows/rnaseqGene/#short-links
#-------------------------------------------------------------------------------------------------------

module load star-rna/2.4.2a  #No need to load since I have already added a SOFT LINK to my bin folder


#----------------------------------------------
#First we generate the genome files
#----------------------------------------------
mkdir GenomeIndex
STAR --runMode genomeGenerate --genomeDir GenomeIndex \
--genomeFastaFiles ../../data/genome/AF70-all-genomic-info/af70_20130423.assembly.fsa \
--sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile ../../data/genome/AF70-all-genomic-info/af70_20130423.gff3 \
--sjdbOverhang 100 --runThreadN 12 --sjdbGTFfeatureExon CDS
#----------------------------------------------------------
#After genome files are generated inside GenomeIndex 
#we can Run STAR for both RNAseq CONDITIONS (PDB and CORN)
#----------------------------------------------------------

#PDB
STAR --genomeDir GenomeIndex \
--readFilesIn ../../data/transcriptome/Aspergillus_AF70_PDB_5D_CAGATCAT_L007_R1_001.fastq \
--runThreadN 12 --outFileNamePrefix aligned-AF70-PDB \
--outSAMattributes All --outFilterMultimapNmax 100 --outFilterMismatchNmax 0 --alignIntronMax 1

--outSAMtype BAM SortedByCoordinate
--quantMode TranscriptomeSAM GeneCounts

#CORN
STAR --genomeDir GenomeIndex \
--readFilesIn ../../data/transcriptome/Aspergillus_AF70_Corn_5D_ACAGTGAT_L007_R1_001.fastq \
--runThreadN 12 --outFileNamePrefix aligned-AF70-CORN \
--outSAMattributes All --outFilterMultimapNmax 100 --outFilterMismatchNmax 0 --alignIntronMax 1
#----------------------------------------------
#Creating BAM files for the Mapped reads
#----------------------------------------------
#PDB
module load samtools
samtools view -bS aligned-AF70-PDBAligned.out.sam \
-o aligned-AF70-PDB.bam
#CORN
samtools view -bS aligned-AF70-CORNAligned.out.sam \
-o aligned-AF70-PDB.bam
 
 #---------------------------------------------------------------------
 #2. RUN DESEQ2 IN R AND RETRIEVE THE MOST SIGNIFICANT READS IN A CSV FILE
 #---------------------------------------------------------------------
 #..............
 #..............
 
 #---------------------------------------------------------------------
 #2.1 Working with the .csv file
 #---------------------------------------------------------------------
 
 #---------------------------------------------------------------------
 #2.1.1 Upregulated 
 #---------------------------------------------------------------------
 cd /panfs/panfs.cluster/scratch/andrese52/Aflavus/results/RNASeq-map
 #Retrieve only UPregulated genes
 awk -F "," '$3>"0" {print $1}' Genes-strongestUPregulation.csv > upregGenes.txt 
 #Eliminate any "" character that could have been created
 sed -i 's/"//g' upregGenes.txt
#Extracting only lines from the gff file containing only the upregulated genes
 cat ../../data/genome/AF70-all-genomic-info/af70_20130423.gff3 | grep -wf upregGenes.txt > gffUpregGenes.txt
 
 #Converting the gffUpregGenes.txt to BED file to extract sequences from the genome fasta file
 #The BED file can be either space delimited or tab delimited, in this case I used tab delimited
 awk -v OFS='\t' '{print $1, $4, $5}' gffUpregGenes.txt > gffUpregGenes.bed
 
 #Extracting the genes that are UPregulated from the AF70 genome
 module load bedtools
 fastaFromBed -fi ../../data/genome/AF70-all-genomic-info/af70_20130423.assembly.fsa -bed gffUpregGenes.bed -fo UpregulatedGenes.fasta
 
 
#------------------------------------------------------------------------
#Blasting the output genes with the Aflatoxin gene clusters AF70 and AF36
#------------------------------------------------------------------------

#AF36 
blastn -task blastn -db /panfs/panfs.cluster/scratch/andrese52/Aflavus/input/AY510455.fasta.db -query UpregulatedGenes.fasta \
-out UpregulatedAF70vsAF36GeneCluster.out -num_threads 8 -word_size 7 \
-outfmt "6 qseqid sseqid score length qlen slen qstart qend sstart send pident nident mismatch positive gapopen gaps qcovhsp sallseqid"

#AF70
blastn -task blastn -db /panfs/panfs.cluster/scratch/andrese52/Aflavus/input/AY510453.fasta.db -query UpregulatedGenes.fasta \
-out UpregulatedAF70vsAF70GeneCluster.out -num_threads 8 -word_size 7 \
-outfmt "6 qseqid sseqid score length qlen slen qstart qend sstart send pident nident mismatch positive gapopen gaps qcovhsp sallseqid"

#How many hits have a 100% identity
awk -F "\t" '$11>=100' UpregulatedAF70vsAF36GeneCluster.out | wc -l  #510
awk -F "\t" '$11>=100' UpregulatedAF70vsAF70GeneCluster.out | wc -l  #527
#.............................................................................................
#If the upregulated genes align perfectly with the Aflatoxin gene cluster, 
#we could say that the upregulated genes have a relationship with the Aflatoxin production
# We can count those genes that align perfectly with the aflatoxin gene cluster with the 
# following command:
#.............................................................................................

awk -F "\t" '$4==$5' UpregulatedAF70vsAF36GeneCluster.out | wc -l #0
awk -F "\t" '$4==$5' UpregulatedAF70vsAF70GeneCluster.out | wc -l #0

 #---------------------------------------------------------------------
 #2.1.2 Downregulated 
 #---------------------------------------------------------------------
 cd /panfs/panfs.cluster/scratch/andrese52/Aflavus/results/RNASeq-map
 #Retrieve only UPregulated genes
 awk -F "," '$3<"0" {print $1}' Genes-strongestUPregulation.csv > DownregGenes.txt 
 #Eliminate any "" character that could have been created
sed -i 's/"//g' DownregGenes.txt
#Eliminate any empty lines (Very important otherwise it retrieves all the lines from the gff lines)
sed -i '/^\s*$/d' DownregGenes.txt

#Extracting only lines from the gff file containing only the upregulated genes
 cat ../../data/genome/AF70-all-genomic-info/af70_20130423.gff3 | grep -wf DownregGenes.txt > gffDownregGenes.txt
 
 #Converting the gffUpregGenes.txt to BED file to extract sequences from the genome fasta file
 #The BED file can be either space delimited or tab delimited, in this case I used tab delimited
 awk -v OFS='\t' '{print $1, $4, $5}' gffDownregGenes.txt > gffDownregGenes.bed
 
 #Extracting the genes that are UPregulated from the AF70 genome
 module load bedtools
 fastaFromBed -fi ../../data/genome/AF70-all-genomic-info/af70_20130423.assembly.fsa -bed gffDownregGenes.bed -fo DownregulatedGenes.fasta
 
 
#------------------------------------------------------------------------
#3. Blasting the output genes with the Aflatoxin gene clusters AF70 and AF36
#------------------------------------------------------------------------

#AF36 
blastn -task blastn -db /panfs/panfs.cluster/scratch/andrese52/Aflavus/input/AY510455.fasta.db -query DownregulatedGenes.fasta \
-out DownregulatedAF70vsAF36GeneCluster.out -num_threads 8 -word_size 7 \
-outfmt "6 qseqid sseqid score length qlen slen qstart qend sstart send pident nident mismatch positive gapopen gaps qcovhsp sallseqid"

#AF70
blastn -task blastn -db /panfs/panfs.cluster/scratch/andrese52/Aflavus/input/AY510453.fasta.db -query DownregulatedGenes.fasta \
-out DownregulatedAF70vsAF70GeneCluster.out -num_threads 8 -word_size 7 \
-outfmt "6 qseqid sseqid score length qlen slen qstart qend sstart send pident nident mismatch positive gapopen gaps qcovhsp sallseqid"

#How many hits have a 100% identity
awk -F "\t" '$11>=100' DownregulatedAF70vsAF36GeneCluster.out | wc -l  #208
awk -F "\t" '$11>=100' DownregulatedAF70vsAF70GeneCluster.out | wc -l  #210

#.............................................................................................
#If the upregulated genes align perfectly with the Aflatoxin gene cluster, 
#we could say that the upregulated genes have a relationship with the Aflatoxin production
# We can count those genes that align perfectly with the aflatoxin gene cluster with the 
# following command:
#.............................................................................................
awk -F "\t" '$4==$5' DownregulatedAF70vsAF36GeneCluster.out >  DownregulatedAF70vsAF36GeneCluster_100ID.out
awk -F "\t" '$4==$5' DownregulatedAF70vsAF36GeneCluster.out | wc -l #14
awk -F "\t" '$4==$5' DownregulatedAF70vsAF70GeneCluster.out > DownregulatedAF70vsAF70GeneCluster_100ID.out
awk -F "\t" '$4==$5' DownregulatedAF70vsAF70GeneCluster.out | wc -l #17
#----------------------------------------------------------------------------
#3.1 Extracting all the sequences that only matched with the AF gene cluster
#----------------------------------------------------------------------------

awk -F "\t" '{print $1}' DownregulatedAF70vsAF70GeneCluster_100ID.out > DownregulatedAF70GeneClusterIDs.txt
module load samtools
xargs samtools faidx DownregulatedGenes.fasta < DownregulatedAF70GeneClusterIDs.txt > DownregulatedAF70GeneCluster.fasta

#------------------------------------------------------------------------
#4. Finding the SNPs to generate e-probes
#------------------------------------------------------------------------
module load mummer
nucmer -maxmatch /panfs/panfs.cluster/scratch/andrese52/Aflavus/input/AY510455.fasta DownregulatedAF70GeneCluster.fasta -p AF70vsAF36Downreg
delta-filter -q AF70vsAF36Downreg.delta > AF70vsAF36Downreg.filter.delta
show-snps -lqT AF70vsAF36Downreg.filter.delta > AF70vsAF36Downreg-SNPs.txt
show-coords AF70vsAF36Downreg.filter.delta -T > AF70vsAF36Downreg.filter.coords

#------------------------------------------------------------------------
#4.1 Retrieving sequences where the SNPs are found as e-probes
#------------------------------------------------------------------------
#...............................................................
#4.1.1  Creating a BED file for the SNPs and e-probes
#...............................................................

#PLEASE NOTE that I did this part using only excel for rapidness, but I could have
#also written a perl or phyton script but since it was gonna be used only once please
#check in the directory RNAseq for E-probe coords-design.xlsx


#------------------------------------------------------------------------
#4.1.1  Retrieving the sequences by using the BED file
#------------------------------------------------------------------------

module load bedtools
fastaFromBed -fi DownregulatedAF70GeneCluster.fasta -bed AF70downreg-eprobe-80.BED -fo AF70-Down-80-eprobes.fasta

#Eliminating redundant e-probe sequences

grep -v '>\|^$'  AF70-Down-80-eprobes.fasta | sort -u > Sequences    # Extract only nucleotide sequences (non-redundant) 
 while read SEQUENCE; do
     grep -w -m 1 -B1 $SEQUENCE AF70-Down-80-eprobes.fasta >> AF70-Down-80-eprobes-uniq.fasta
 done < Sequences

#HERE YOU RUN EDNA SCRIPT WITH THE DESIGNED E-PROBES AND CAN COUNT HOW MANY READS ALIGNED TO EACH PROBE BY RUNNING:
#THIS CODE SHOULD RUN IN THE OUTPUT FOLDER OF EDNA OR ednaout
awk '{print $2}' AF70-Corn-flavus-AF70-toxin-80-parsed.txt-raw | uniq --count | sort
  
  