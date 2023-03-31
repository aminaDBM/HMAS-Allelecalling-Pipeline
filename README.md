# HMAS-Allelecalling-Pipeline

* [Introduction](#introduction)
* [Pipeline Workflow](#pipeline-workflow)

### Introduction
HMAS-Allelecalling pipeline is designed to identify alleles in metagenomics samples using highly-multiplexed amplicon sequencing data. This pipeline does reference based BLAST search for allele calling. It acts as a denoising step for being able to identify different variants of Salmonella with higher accuracy in the downstream analysis.

![Screenshot (27)](https://user-images.githubusercontent.com/93733968/228165014-f6c069b0-66eb-4478-832d-149f284869f7.png)


### Pipeline Workflow
Allele calling pipeline consists of four integral steps defined as follows:

![Screenshot (32)](https://user-images.githubusercontent.com/93733968/228159991-b35b6d77-dac5-434d-b569-1eac641194f2.png)

#### Step1:

Run BLAST on input "fasta" files.
* INPUT: fasta files, InsilicoPrimers database
* OUTPUT: tsv files

#### Step2:

Run Rscript "Known.alleles-pipelinepart1.R" to apply filters on the first BLAST results and identify "Known Alleles".
* INPUT: fasta files, tsv files
* OUTPUT: KnownAlleles.fasta, KnownAlleleStatistics.csv, KnownAlleles.consolidatedStats.csv, UnknownAlleles.fasta

#### Step3:

Run BLAST on "Unknown Alleles" fasta files
* INPUT: UnknownAlleles fasta files, Centroid Reference database
* OUTPUT: UnknownAlleles.tsv files

#### Step4:

Run Rscript "Novel.alleles-pipelinepart2.R" to apply filters on second BLAST results and identify "Novel Alleles".
* INPUT: UnknownAlleles fasta files, UnknownAlleles.tsv files
* OUTPUT: NovelAlleles.fasta, NovelAlleleStatistics.csv, NovelAlleles.consolidatedStats.csv, NotfoundAlleles.fasta
