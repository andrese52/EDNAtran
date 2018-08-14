#!/usr/bin/env bash

#Download all the sequences required

perl scripts/get_genbank_batch.pl $1
perl scripts/genbank2fasta.pl sequences.gbk
