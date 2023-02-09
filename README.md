## Centroid_Reference Generation
----------------------------------

This code has been accredited to:

BacTag - pipeline for fast and accurate genes and alleles typing (https://git.lumc.nl/l.khachatryan/BacTag.git)

* [Preparations](#preparations)
* [Installation](#installation)
* [Centroid Reference Generation](#centroid reference generation)


## Preparations

-----------------

Creating a conda env using yml file

1. Create the environment from the allelecalling.yml file:

      conda env create -f allelecalling.yml
2. Activate the new environment: 
      conda activate allelecalling
3. Verify that the new environment was installed correctly:
     conda env list

## Installation
-----------------
From source:
Centroid_Reference Generation 
   git clone https://github.com/aminaDBM/Centroid_Reference.git
Fastools
  git clone https://git.lumc.nl/j.f.j.laros/fastools.git
  cd fastools
  pip install .
Sim-reads
  git clone https://git.lumc.nl/j.f.j.laros/sim-reads.git
  cd sim-reads
  pip install .

##Database setup
-------------------
###Database organisation requirements
*  The database should be placed in `data` folder.
*  It should consist of allele sequences in
  [FASTA](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHel)
  format saved in a file with `.fa` file extension. 
*  The name of each allele sequence should contain the gene name and allele
  name, separated by `-`. 
*  Alleles of the same gene should be placed together to the folder with gene
  name.

