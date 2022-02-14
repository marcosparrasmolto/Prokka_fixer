#!/bin/bash

DIR="/Output/"

#for i in ./Acces_genomes/*fa;do
if [ -d "$DIR" ]; then
  # Take action if $DIR exists. #
  mkdir Output
fi


cd ./Seqs/
for i in *fa;do
#perl ../Juntar_fasta.pl $i
prodigal -i "$i" -d aux_file -a aux_file_aa -q
cat aux_file >> ../aux_file_prodigal_global
cat aux_file_aa >> ../aux_file_prodigal_global_aa
blastn -max_hsps 1 -query aux_file -outfmt 6 -db ../Database/ETEC_DB -out salida_blast
cat salida_blast >> ../blast_output_file
rm -rf salida_blast
prokka "$i" --outdir ../Output/"$i"_prokka_out --prefix "$i" -rfam --locustag ETEC --force
done

#--norrna --notrna

cd ../

grep ">" aux_file_prodigal_global > aux_file_prodigal_global_headers.txt
