#!/bin/bash


############ BUILDING DATABASES ###########
#Description: Allelecalling pipeline rrequires two databases to be constructed for running two BLAST runs in sorting sequences
#as Known alleles and Novel alleles.

#Database for first BLAST run (BLAST1_DB)
echo "Provide the fasta file for building database for first BLAST run: "
read BLAST1fasta
echo "Title for database: "
read BLAST1title
echo "Path to save database for first BLAST run: "
read BLAST1dir

makeblastdb -in $BLAST1fasta -title $BLAST1title -dbtype nucl -out $BLAST1dir/$BLAST1title


#Database of Centroid Reference for second BLAST run (BLAST2_DB)
echo "Provide the fasta file for building centroid reference database for second BLAST run: "
read BLAST2fasta
echo "Title for database: "
read BLAST2title
echo "Path to save centroid reference database for second BLAST run: "
read BLAST2dir

makeblastdb -in $BLAST2fasta -title $BLAST2title -dbtype nucl -out $BLAST2dir/$BLAST2title

###########################################

