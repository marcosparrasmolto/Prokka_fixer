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
prodigal -i "$i" -d salida_prueba -a salida_prueba_aa -q
cat salida_prueba >> ../salida_prodigal_global
cat salida_prueba_aa >> ../salida_prodigal_global_aa
blastn -max_hsps 1 -query salida_prueba -outfmt 6 -db ../Database/ETEC_DB -out salida_blast
cat salida_blast >> ../salida_global_blast
prokka "$i" --outdir ../Output/"$i"_prokka_out --prefix "$i" -rfam --force
done

#--norrna --notrna

cd ../

grep ">" salida_prodigal_global > salida_prodigal_global_headers.txt
