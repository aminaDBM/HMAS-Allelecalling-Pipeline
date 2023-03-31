#!/bin/bash

echo "input fasta file"
read Inputfasta

echo "Gene list text file"
read GeneList

demuxbyname.sh in=$Inputfasta  out=%.fasta names=$GeneList substringmode

for x in ./*.fasta; do
  mkdir "${x%.*}" && mv "$x" "${x%.*}"
done

for f in */; do ( cd "$f"; awk -F "|" '/^>/ {close(F) ; F = substr($1,2,length($1)-1)".fasta"} {print >> F}' ${f%/}.fasta ); done

for f in */; do (cd ${f%/}; rm ${f%/}.fasta); done
