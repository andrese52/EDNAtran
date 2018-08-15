#!/usr/bin/env bash

#Download all the sequences required

get_genbank_batch.pl $1
genbank2fasta.pl sequences.gbk
