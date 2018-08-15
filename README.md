Running EDNAtran with *Aspergillus flavus* as a model system.
===
Contents
[TOC]

This repository contains all the [scripts](/scripts) that allowed to obtain all the [231 e-probes](/results/eprobes/AF70-80.fasta) associated with aflatoxin production in AF70.
There are multiple steps associated with the design of e-probes. In this case, all the methodology is explained in our manuscript.

# E-probe coordinates
Retrieving coordinates from different aflatoxin gene clusters. First we need to retrieve all the aflatoxin gene clusters available on NCBI.
First move to the `results` folder and run the following command

```bash
cd results
bash ../scripts/run_circos.sh ../gene-cluster-IDs.txt
```
