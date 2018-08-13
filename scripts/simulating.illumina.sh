#!/bin/bash


#------------------------------------------------------------------------------#
# Obtain an updated and detailed list of the assemblies directories in the NCBI FTP site
#See this ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt and this 
#https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#asmsumfiles
#------------------------------------------------------------------------------#
#================Aflavus=============
cd flavus/ftp.ncbi.genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
awk -F'\t' '$8~"Aspergillus flavus" && $14 == "Full" && $13 =="Major"' assembly_summary_genbank.txt >  ftp-flavus-assemblies.inf.txt
awk -F'\t' '{print $20}' ftp-flavus-assemblies.inf.txt > ftp-flavus-assemblies.txt
#Downloading the ftp folders in parellel
cat ftp-flavus-assemblies.txt | parallel --gnu "wget -m {}"
#You will get the folder directories organized so then we have to extract each fna and gff files to another folder
find ftp.ncbi.nlm.nih.gov/ -name '*.fna*' -exec cp {} $HOME/projects/flavus/ \;
find ftp.ncbi.nlm.nih.gov/ -name '*.gff*' -exec cp {} $HOME/projects/flavus/ \;
cd $HOME/projects/flavus
rm *cds*
rm *rna*
gunzip *.gz
#================Corn=============
awk -F'\t' '$8~"Zea mays" && $14 == "Full" && $13 =="Major"' assembly_summary_genbank.txt >  ftp-zeamays-assemblies.inf.txt
awk -F'\t' '{print $20}' ftp-zeamays-assemblies.inf.txt > ftp-zeamays-assemblies.txt
#Downloading the ftp folders in parellel
cat ftp-zeamays-assemblies.txt | parallel --gnu "wget -m {}"
#You will get the folder directories organized so then we have to extract each fna and gff files to another folder
find ftp.ncbi.nlm.nih.gov/ -name '*.fna*' -exec cp {} $HOME/projects/flavus/genomes \;
find ftp.ncbi.nlm.nih.gov/ -name '*.gff*' -exec cp {} $HOME/projects/flavus/genomes \;
cd $HOME/projects/flavus/genomes
rm *cds*
rm *rna*
gunzip *.gz




#============Multifasta to Single fasta========================

#Converting to single fasta
ls *.fna > list.fasta.txt
while read F  ; do
        awk '
    /^>/ { 
        # print the first header
        if (c++ == 0) {print; print ""} 
        next
    } 
    /^$/ {next} 
    {printf "%s", $0} 
    END {print ""}
' $F > $F-singlefa.fasta
done < list.fasta.txt
#--------------------------------

head -n 1 *.fasta > headers.txt
#Adding the fasta files to databases
MetaSim cmd --add-files *.fasta

#Simulating 454
MetaSim cmd -4 --454-paired-read-length 20 --454-mate-probability 0 --454-logn-mean 0.23 --454-cycles 200 --454-logn-std 0.15 --454-nosqrt -r 5000 -d /home/andres_espindolac/projects/flavus/simulations/flavus-consortium.fasta -z headers.mprf

#Simulating Illumina from here https://github.com/afbase/Pyllumina/blob/master/make_reads_Illumina.sh
MetaSim cmd -r 200000 -m \
-g Illumina-100.mconf -2 Illumina-100.mconf \
--empirical-pe-probability 100 --clones-mean 350 --clones-param2 800 \
-d /home/andres_espindolac/projects/flavus/simulations/flavus-consortium.fasta -z headers.mprf > metasim.log
