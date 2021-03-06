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
  show-snps -TH $i.filter.delta > $i.snps
  coords-nucmer-parser-circos.pl $i.coords $i
  makeblastdb -in $i -dbtype nucl
  blastn -db $i -query eprobes/AF70-80.fasta -outfmt 6 -out $i-vs-AF70-80.out
  awk '{print $2,$9,$10}' $i-vs-AF70-80.out > $i-PROBE-circos-highlight.txt
  awk '{print $2,$9,$10,$1}' $i-vs-AF70-80.out > $i-PROBE-circos-bands.txt
done
cat sequences.fasta.files/*PROBE* > highlight-PROBE-circos.txt
cat sequences.fasta.files/*.snps > all-snps.txt
cat sequences.fasta.files/*circos-bands.txt > all-circos-band.txt
ln -s $(pwd)/highlight-PROBE-circos.txt $(pwd)/circos-af/data/
ln -s $(pwd)/all-circos-band.txt $(pwd)/circos-af/data/
