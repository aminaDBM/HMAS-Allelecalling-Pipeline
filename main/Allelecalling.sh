#!/bin/bash

#Taking path to data from user
echo "Enter the path for data: "
read path



ls $path*_extractedAmplicons.fasta | parallel -a - blastn -query {} -db /scicomp/home-pure/spu0/19isolatesallprimersdata/sal19_extractedAmplicon -perc_identity 100 -outfmt \"6 qseqid sseqid pident qcovs qlen slen length \" -out {.}.tsv


Rscript Known.alleles-pipelinepart1.R $path*_extractedAmplicons.tsv $path*_extractedAmplicons.fasta |& tee -a Known.alleles-pipelinepart1.log



#ls *.unknown_Alleles.fasta | parallel -a - blastn -query {} -db /scicomp/home-pure/spu0/266Pacbio.IlluminaData/sal19.Centroid_RefDB -outfmt \"6 qseqid sseqid pident qcovs gaps btop qlen slen length \" -out {.}.tsv


#Rscript Unknown.alleles-pipelinepart2.R *.unknown_Alleles.tsv .unknown_Alleles.fasta |& tee -a Unknown.alleles-pipelinepart2.log


