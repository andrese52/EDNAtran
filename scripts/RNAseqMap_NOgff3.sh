#!/bin/bash


###THIS OPTION DOES NOT USE GFF3 FILES


##Testing without gff3 files
STAR --runMode genomeGenerate --genomeDir GenomeIndex \
--genomeFastaFiles ../../data/genome/AF70-all-genomic-info/af70_20130423.assembly.fsa \
--runThreadN 12

#After genome files are generated inside GenomeIndex we can Run STAR for both RNAseq CONDITIONS (PDB and CORN)
#PDB
STAR --genomeDir GenomeIndex \
--readFilesIn ../../data/transcriptome/Aspergillus_AF70_PDB_5D_CAGATCAT_L007_R1_001.fastq \
--runThreadN 12 --outFileNamePrefix aligned-AF70-PDB 

#CORN
STAR --genomeDir GenomeIndex \
--readFilesIn ../../data/transcriptome/Aspergillus_AF70_Corn_5D_ACAGTGAT_L007_R1_001.fastq \
--runThreadN 12 --outFileNamePrefix aligned-AF70-CORN

#Creating BAM files for the Mapped reads
#PDB
module load samtools
samtools view -bS aligned-AF70-PDBAligned.out.sam \
-o aligned-AF70-PDB.bam
#CORN
samtools view -bS aligned-AF70-CORNAligned.out.sam \
-o aligned-AF70-CORN.bam

#From Here and On continue with Bioconductor and R

 