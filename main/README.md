## HMAS-Allelecalling Pipeline
![Screenshot (26)](https://user-images.githubusercontent.com/93733968/228119258-6d765bd0-3294-4b54-859e-b51e38083eef.png)

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

Run Rscript "Unknown.alleles-pipelinepart2.R" to apply filters on second BLAST results and identify "Novel Alleles".
* INPUT: UnknownAlleles fasta files, UnknownAlleles.tsv files
* OUTPUT: NovelAlleles.fasta, NovelAlleleStatistics.csv, NovelAlleles.consolidatedStats.csv, NotfoundAlleles.fasta

                                                             

