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
  show-coords $i.filter.delta -qcl > $i.coords
  coords-nucmer-parser-circos.pl $i.coords $i
  makeblastdb -in sequences.fasta.files/$i -dbtype nucl
  blastn -db $i -query eprobes/AF70-80.fasta -outfmt 6 -out $i-vs-AF70-80.out

done
