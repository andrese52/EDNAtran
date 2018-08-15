#!/usr/bin/env bash

#Download all the sequences required

perl get_genbank_batch.pl $1
perl genbank2fasta.pl sequences.gbk
