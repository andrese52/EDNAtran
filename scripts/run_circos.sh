#!/usr/bin/env bash

#Download all the sequences required

get_genbank_batch.pl $1
genbank2fasta.pl sequences.gbk
separatemultifasta.pl sequences.fasta
# RUNNING NUCMER
#------------------------------------------------------------------------------#

#nucmer -maxmatch ref.fasta qry.fasta -p output-name
#delta-filter -q output-name.delta > output-name.filter.delta
#show-coords output-name.filter.delta -T > output-name.coords
