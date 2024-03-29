#!/bin/bash

######################INPUTS#################

while getopts i:f:s: flag
do
    case "${flag}" in
        i) path=${OPTARG};;
        f) BLAST1_DB=${OPTARG};;
        s) BLAST2_DB=${OPTARG};;
    esac
done

#USAGE: bash Allelecalling.sh -i /path/to/input/data -f /name/and/path/to/first/BLAST/database/BLAST1_DB -s /name/and/path/to/second/BLAST/database/BLAST2_DB 

#Description: This section of the script takes in path as input and stores them as variables:

#1. $path:  Path to input data files
#2. $BLAST1_DB: Path and title to the database used to run first BLAST for preliminary sorting of input sequences as "Known Alleles" and "Unknown Alleles"
#3. $BLAST2_DB: Path and title to the Centroid reference database for running second BLAST to sort "Unknown Alleles" into "Novel Alleles" and "Alleles Not Found"




################# PIPELINE STEPS ###################

#Description: This section contains FOUR execution commands for the Allelecalling pipeline. The four processes can be defined as follows:

#STEP1:  

#Run BLAST on input "fasta" files.
#INPUT: fasta files, InsilicoPrimers database
#OUTPUT: tsv files

#STEP2: 

#Run Rscript "Known.alleles-pipelinepart1.R" to apply filters on the first BLAST results and identify "Known Alleles".
#INPUT: fasta files, tsv files
#OUTPUT: KnownAlleles.fasta, KnownAlleleStatistics.csv, KnownAlleles.consolidatedStats.csv, UnknownAlleles.fasta

#STEP3: 

#Run BLAST on "Unknown Alleles" fasta files
#INPUT: UnknownAlleles fasta files, Centroid Reference database
#OUTPUT: UnknownAlleles.tsv files

#STEP4: 

#Run Rscript "Novel.alleles-pipelinepart2.R" to apply filters on second BLAST results and identify "Novel Alleles".
#INPUT: UnknownAlleles fasta files, UnknownAlleles.tsv files
#OUTPUT: NovelAlleles.fasta, NovelAlleleStatistics.csv, NovelAlleles.consolidatedStats.csv, NotfoundAlleles.fasta

							     ########################################################

############### EXECUTION COMMANDS #####################


#1. Run BLAST on input "fasta" files. (INPUT: fasta files, InsilicoPrimers database) (OUTPUT: tsv files)

ls $path/*_extractedAmplicons.fasta | parallel -a - blastn -query {} -db $BLAST1_DB  -perc_identity 100 -outfmt \"6 qseqid sseqid pident qcovs qlen slen length \" -out {.}.tsv



#2. Run Rscript "Known.alleles-pipelinepart1.R" to apply filters on the first BLAST results and identify "Known Alleles". (INPUT: fasta files, tsv files) (OUTPUT: KnownAlleles.fasta, KnownAlleleStatistics.csv, KnownAlleles.consolidatedStats.csv, UnknownAlleles.fasta)

export Path_to_Input=$path
echo -e 'source("Known.alleles-pipelinepart1.R")  \n q()' | R --no-save --slave


# Creating directory for Known Alleles Output data
OutputDir="$path/KnownAlleles_Output"
#create output directory
mkdir -p  $OutputDir
mv $path/*.known_Alleles.fasta $path/*-KnownAlleleStatistics.csv $path/KnownAlleles.consolidatedStats.csv $OutputDir

cat $OutputDir/*.known_Alleles.fasta* > $path/AllKnownAlleles.fasta

#Creating directory for Unknown Alleles.fasta files (subject to further processing)
Unknown_AllelesInput="$path/UnknownAlleles_Data"
#create output directory
mkdir -p  $Unknown_AllelesInput
mv $path/*unknown_Alleles.fasta  $Unknown_AllelesInput



#3. Run BLAST on "Unknown Alleles" fasta files (INPUT: UnknownAlleles fasta files, Centroid Reference database) (OUTPUT: UnknownAlleles.tsv files)

ls $Unknown_AllelesInput/*.unknown_Alleles.fasta | parallel -a - blastn -query {} -db $BLAST2_DB -outfmt \"6 qseqid sseqid pident qcovs gaps btop qlen slen length \" -out {.}.tsv


#4. Run Rscript "Novel.alleles-pipelinepart2.R" to apply filters on second BLAST results and identify "Novel Alleles". (INPUT: UnknownAlleles fasta files, UnknownAlleles.tsv files) (OUTPUT: NovelAlleles.fasta, NovelAlleleStatistics.csv, NovelAlleles.consolidatedStats.csv, NotfoundAlleles.fasta) 

export Path_to_Input=$Unknown_AllelesInput
echo -e 'source("Novel.alleles-pipelinepart2.R")  \n q()' | R --no-save --slave

#awk 'FNR==1{print ">" FILENAME; next}1' *.Novel_Alleles.fasta > AllNovelAlleles.fasta



#Creating directory for Novel Alleles data
Novel_AllelesOutput="$path/NovelAlleles_Output"
#create output directory
mkdir -p  $Novel_AllelesOutput
mv $Unknown_AllelesInput/*.Novel_Alleles.fasta $Unknown_AllelesInput/*.NovelAlleleStatistics.csv $Unknown_AllelesInput/NovelAlleles.consolidatedStats.csv $Novel_AllelesOutput

#Creating directory for Notfound Alleles data
Notfound_AllelesOutput="$path/NotfoundAlleles_Output"
#create output directory
mkdir -p  $Notfound_AllelesOutput
mv $Unknown_AllelesInput/*.Notfound_Alleles.fasta $Notfound_AllelesOutput

cat $Novel_AllelesOutput/*.Novel_Alleles.fasta* > $path/AllNovelAlleles.fasta


#Output file for strainRecon by concatenating all known alleles and all novel alleles

cat $path/AllKnownAlleles.fasta $path/AllNovelAlleles.fasta > $path/InputforStrainRecon.fasta





