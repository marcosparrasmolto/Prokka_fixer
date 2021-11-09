
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **README**

## **R code for: Parras-Moltó & , 2021**

This document describes the use of Prokka\_anotation\_fixer.R. The
analysis of ETEC E.coli species could be a problem for certain proteins
due to their similarity to not related pathogenicity proteins form in
Porkka database, so we could get a biased result when using default the
defaul one. A blast using a custom ETEC database could help us to fix
some ETEC pathogenicity related genes and modify the Prokka output with
new, custom, results.

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

3- `./Databse/ETEC_DB`. Custom ETEC proteins database.

**Usage**

Multifasta assembly/genomes files should be stored in a folder called
*Seqs* in the same folder this script is being run.

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
## 2  ERS055657.7114_1_10.27_32
## 3   ERS055657.7114_1_10.59_4
## 4   ERS055657.7114_1_10.59_5
## 5   ERS055657.7114_1_10.59_6
## 6   ERS055657.7114_1_10.59_7
##                                                                            sseqid
## 1                      NZ_AFAH00000000.2_fim41a-G_F41a_Fim41a-G_fimbrial_subunit_
## 2 CBJ02741.1_yghJ_YghJ_secreted_lipoprotein,_mucine_binding_metalloprotease_YghJ_
## 3                                              LR883052_lngR_CS21_LngR_regulator_
## 4                                              LR883052_lngS_CS21_LngS_regulator_
## 5                                       LR883052_lngT_CS21_LngT_unknown_function_
## 6                                      LR883052_lngX2_CS21_LngX2_unknown_function
##    pident length mismatch gapopen qstart qend sstart send    evalue bitscore
## 1 100.000     53        0       0     89  141    354  302  4.41e-23       99
## 2  86.874   4594      541      48      1 4566      1 4560  0.00e+00     5086
## 3 100.000    306        0       0      1  306      1  306 2.33e-163      566
## 4 100.000    897        0       0      1  897      1  897  0.00e+00     1657
## 5 100.000    444        0       0      1  444      1  444  0.00e+00      821
## 6 100.000    222        0       0      1  222      1  222 8.24e-117      411
##                      ACC  Protein                          GCA
## 1 ERS055657.7114_1_10.14 fim41a-G 7114_1_10.contigs.fa_OUT.txt
## 2 ERS055657.7114_1_10.27     YghJ 7114_1_10.contigs.fa_OUT.txt
## 3 ERS055657.7114_1_10.59     CS21 7114_1_10.contigs.fa_OUT.txt
## 4 ERS055657.7114_1_10.59     CS21 7114_1_10.contigs.fa_OUT.txt
## 5 ERS055657.7114_1_10.59     CS21 7114_1_10.contigs.fa_OUT.txt
## 6 ERS055657.7114_1_10.59     CS21 7114_1_10.contigs.fa_OUT.txt
```

2- `Composition_by_assembly.txt`. Will save the information about the
proteins included for each genome.

``` r
head(composition_genome)
##                             V1
## 1 7114_1_10.contigs.fa_OUT.txt
##                                                                         V2
## 1 fim41a-G, YghJ, eatA, etpA, etpB, eltA, etpC, CFA/I, STh, CS21, CS1, CS3
```

3- `./Output/*gff`. Modified Prokka .gff files will be stored in Output
folder with the string “new\_gff\_file” at the end of file name, for
each analyzed file.
