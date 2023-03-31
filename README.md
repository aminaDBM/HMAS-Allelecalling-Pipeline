# HMAS-Allelecalling-Pipeline

* [Introduction](#introduction)
* [Pipeline Workflow](#pipeline-workflow)

### Introduction
HMAS-Allelecalling pipeline is designed to identify alleles in metagenomics samples using highly-multiplexed amplicon sequencing data. This pipeline does reference based BLAST search for allele calling. It acts as a denoising step for being able to identify different variants of Salmonella with higher accuracy in the downstream analysis.

![Screenshot (39)](https://user-images.githubusercontent.com/93733968/229218647-c54b95c6-65bc-4d5f-912a-f4f70a1defac.png)

### Pipeline Workflow
Allele calling pipeline consists of four integral steps defined as follows:

![Screenshot (40)](https://user-images.githubusercontent.com/93733968/229219275-db18fcdc-ed7f-4f93-82d6-5e1b889d221d.png)
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
