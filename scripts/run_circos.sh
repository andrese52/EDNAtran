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
  delta-filter -q $i.nucmer.delta > $i.filter.delta
  show-coords $i.filter.delta -Trcl > $i.coords
  coords-nucmer-parser-circos.pl $i.coords $i
done
