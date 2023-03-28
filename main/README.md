## HMAS-Allelecalling Pipeline
![Screenshot (27)](https://user-images.githubusercontent.com/93733968/228119977-c15fa1ab-7025-4b26-b957-daeef87e0622.png)
### Pipeline Steps 

Allelecalling pipeline consists of four integral steps defined as follows:

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

### Running the pipeline

The pipeline comprises of Rscripts and a bash script kept in the same path. The bash script requires the user to specigy paths to input data, Insilico primers database and centroid reference database.

  bash Allelecalling.sh
  
