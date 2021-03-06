
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **README**

## **R code for: Parras-Moltó & , 2021**

This document describes the use of Prokka\_annotation\_fixer.R. Annotation of enterotoxigenic Escherichia coli (ETEC) genomes can be problematic. Certain ETEC virulence proteins are similar to proteins which are not linked to pathogenicity. A blast using a custom ETEC database with ETEC virulence factors with this script will help generate annotated genomes by modifying the Prokka output with new, custom, annotations.

    R version 4.1.0 (2021-05-18)

### **Dependencies**

#### **Conda**

A Conda environment is required to run the script.

Miniconda3 installation guide:

``` r
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
```

Environment installation:

``` r
conda env create -f prokka_env.yml
```

**Packages**

· `stringr` v1.4.0. A consistent, simple and easy to use set of wrappers
around the fantastic ‘stringi’ package. All function and argument names
(and positions) are consistent, all functions deal with “NA”’s and zero
length vectors in the same way, and the output from one function is easy
to feed into the input of another.

**Files located in the same folder as the script**

1- `./bin/blaster.sh`. Script needed to run prodigal, blast and prokka.
Results are generated at this point..

2- `./bin/prokka_parser.pl`. Script that parses all the information from
the original Prokka file and change the annotation for a new one
produced from blast against the custom database.

3- `./Database/ETEC_DB`. Custom ETEC proteins database.

4- `./Database/ETEC_CFs`. Criteria for ETEC classification.

How to create the database:

· Headers from the fasta file should be named like:

``` r
>id-?gene?product?-
```

· A database is created executing the folling command line:

``` r
makeblastdb -in database.fasta -dbtype nucl -out ETEC_DB
```

**Usage**

Multifasta assembly/genomes files should be stored in a folder called
*Seqs* in the same folder this script is being run. Files should be as
\*.fa format.

First of all you will need to load the prokka environment in the bash
terminal

``` r
conda activate prokka_env
```

After that, execute the script:

``` r
Rscript Prokka_anotation_fixer.R
```

***Output results are found in three different files:***

1- `All_results_hit.txt`. Will save the best hit for each ORF based on
the most abundant one.

``` r
head(best_blast)
##                       qseqid
## 1 ERS055657.7114_1_10.14_111
## 2   ERS055657.7114_1_10.59_4
## 3   ERS055657.7114_1_10.59_5
## 4   ERS055657.7114_1_10.59_6
## 5   ERS055657.7114_1_10.59_7
## 6   ERS055657.7114_1_10.59_8
##                                                        sseqid pident length
## 1 NZ_AFAH00000000.2-?fim41a-G?F41a-Fim41a-G-fimbrial-subunit?    100     53
## 2                         LR883052-?lngR?CS21-LngR-regulator?    100    306
## 3                         LR883052-?lngS?CS21-LngS-regulator?    100    897
## 4                  LR883052-?lngT?CS21-LngT-unknown-function?    100    444
## 5                 LR883052-?lngX2?CS21-LngX2-unknown-function    100    222
## 6                     LR883052-?lngA?CS21-LngA-major-subunit?    100    711
##   mismatch gapopen qstart qend sstart send    evalue bitscore
## 1        0       0     89  141    354  302  4.06e-23       99
## 2        0       0      1  306      1  306 2.15e-163      566
## 3        0       0      1  897      1  897  0.00e+00     1657
## 4        0       0      1  444      1  444  0.00e+00      821
## 5        0       0      1  222      1  222 7.59e-117      411
## 6        0       0      1  711      1  711  0.00e+00     1314
##                      ACC Protein                  GCA
## 1 ERS055657.7114_1_10.14    F41a 7114_1_10.contigs.fa
## 2 ERS055657.7114_1_10.59    CS21 7114_1_10.contigs.fa
## 3 ERS055657.7114_1_10.59    CS21 7114_1_10.contigs.fa
## 4 ERS055657.7114_1_10.59    CS21 7114_1_10.contigs.fa
## 5 ERS055657.7114_1_10.59    CS21 7114_1_10.contigs.fa
## 6 ERS055657.7114_1_10.59    CS21 7114_1_10.contigs.fa
```

2- `Composition_by_assembly.txt`. Will save the information about the
proteins included for each genome.

``` r
head(composition_genome)
##                     V1                       V2
## 1 7114_1_10.contigs.fa F41a, CS1, CS3, LTh, STh
```

3- `./Output/*gff`. Modified Prokka .gff files will be stored in Output
folder with the string “new\_gff\_file” at the end of file name, for
each analyzed file.
