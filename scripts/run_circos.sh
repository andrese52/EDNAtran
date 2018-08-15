#!/usr/bin/env bash

#Download all the sequences required

get_genbank_batch.pl $1
genbank2fasta.pl sequences.gbk
separatemultifasta.pl sequences.fasta
# RUNNING NUCMER
#------------------------------------------------------------------------------#
for i in sequences.fasta.files/*.fasta
do
  echo $i
  nucmer -maxmatch sequences.fasta.files/AY510453.fasta $i -p $i.nucmer
done

#delta-filter -q output-name.delta > output-name.filter.delta
#show-coords output-name.filter.delta -T > output-name.coords
