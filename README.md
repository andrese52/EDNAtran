Running EDNAtran with *Aspergillus flavus* as a model system.
===
Contents


This repository contains part of the [scripts](/scripts) that allowed to validate all the [231 e-probes](/results/eprobes/AF70-80.fasta) associated with aflatoxin production in AF70.
There are multiple steps associated with the design of e-probes. In this case, all the methodology is explained in our manuscript.

To start replicating the results of our manuscript first clone the repository as follows:

```bash
git clone https://github.com/andrese52/EDNAtran.git
cd EDNAtran
```
# Dependencies
- [Mummer 3.9.4 alpha](https://github.com/mummer4/mummer/releases/tag/v3.9.4alpha): Installed following [these instructions](https://github.com/mummer4/mummer/blob/master/INSTALL.md)
- Bioperl Search module
    ```bash=1
    cpan Bio::SearchIO
    ```
- [Circos](http://circos.ca/software/installation/)
    ```bash=1
    wget http://circos.ca/distribution/circos-0.69-6.tgz
    tar xvfz circos-0.69-6.tgz
    cd circos-0.69-6
    ```
    - Some perl modules are required by circos and can be installed as follows:

    ```bash=1
    sudo apt-get install libgd-perl
    cpan Config::General
    cpan GD::Polyline
    cpan Math::Bezier
    cpan Math::Round
    cpan Readonly
    cpan Params::Validate
    cpan Math::VecStat
    cpan Statistics::Basic
    cpan Regexp::Common
    cpan Text::Format
    cpan Set::IntSpan
    ```
   However, it depends on your OS, in other cases you might need to install more or less modules. This was our specific case. We are using Ubuntu 18.04.
    - You can also download the examples folder to practice [from here](http://circos.ca/distribution/circos-course-2017.tgz)
    Then add the bin folder to your permanent `$PATH`
- gnu-plot
    ```bash=
     sudo apt install gnuplot-x11
    ```
-

# Setting up environmental variables
We need to tell our computer where to find the scripts. Therefore we must run this code once we have cloned the repository:

```bash
export PATH=$PATH:$(pwd)/scripts
```
# E-probe coordinates
Retrieving coordinates from different aflatoxin gene clusters. First we need to retrieve all the aflatoxin gene clusters available on NCBI.
First move to the `results` folder and run the following command

```bash=
cd results
bash ../scripts/run_circos.sh ../gene-cluster-IDs.txt
```
