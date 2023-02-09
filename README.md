## Centroid_Reference Generation
----------------------------------

This code has been accredited to:

BacTag - pipeline for fast and accurate genes and alleles typing (https://git.lumc.nl/l.khachatryan/BacTag.git)

* [Preparations](#preparations)
* [Installation](#installation)
* [Database Setup](#databasesetup)
* [Centroid Reference Generation](#centroidreferencegeneration)


## Preparations
Creating a conda environment using `.yml` file

* Create the environment from the allelecalling.yml file:

      conda env create -f allelecalling.yml
      
* Activate the new environment: 
      
      conda activate allelecalling
      
* Verify that the new environment was installed correctly:
      
      conda env list

## Installation
From source:
* Centroid_Reference Generation 

      git clone https://github.com/aminaDBM/Centroid_Reference.git
      
* Fastools

      git clone https://git.lumc.nl/j.f.j.laros/fastools.git
      cd fastools
      pip install .
      
* Sim-reads

      git clone https://git.lumc.nl/j.f.j.laros/sim-reads.git
      cd sim-reads
      pip install .


## Database Setup
#### Database organization requirements
*  The database should be placed in `data` folder.
*  It should consist of allele sequences in
  [FASTA](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHel)
  format saved in a file with `.fa` file extension. 
*  The name of each allele sequence should contain the gene name and allele
  name, separated by `-`. 
*  Alleles of the same gene should be placed together to the folder with gene
  name.
*  Gene and allele names should consist only of letters and numbers.
#### Example
You have a database (`MLST`), that includes allele sequences of 2 genes
(`GENE1` and `GENE2`), 3 alleles per gene. Thus you have 6 sequences in your
database in total.

![Screenshot (21)](https://user-images.githubusercontent.com/93733968/217869081-bb510c10-7b69-4d37-ad5a-bb2af475330d.png)

The example of correctly organised database (test_data) is already placed to the
`data` folder.

### Configuration file

It should be placed to the `analysis` directory. There are required and
optional fields.

#### Required variables
* `path` - path to the raw database.
#### Optional variables
* `flanks` - file with flanking sequences. If nothing is specified, poli-N
  sequences of length 50 will be used as flanks from both sides.

## Centroid Reference Generation
Before starting the calculations make sure that you have:
1. Raw database In the [proper format](#database-setup) stored in
  the `test_data` folder.
2. [Config file](#configuration-file) with all variables
  provided.

Database processing has to be initiated in the `analysis` folder, the
processed database is going to be stored in the same folder. To initiate a
database processing, type:

    bash build_database.bh

This will launch the 3 step procedure.

#### Step 1
First, the pipeline checks whether the current processed database (if there
is one in this folder) corresponds to the one that user want to process.
Sequences in the processed database that are not present in the raw database
will be deleted, as well as all of the files dependent on them. If there are
sequences in the raw database, that are not present in the processed
database, the processing will be redone for the gene containing new
sequence.

#### Step 2
In the second step, each sequence from the raw database will be merged with
flanking sequences. Flanking sequences can be provided by the user (see
[Flanking regions](#flanking-regions)) or can be the default
fifty-nucleotide long poly-N sequence. After that, the pipeline calculates a
pairwise alignment distance matrix for each gene and provide the quality
reports. The distance matrix will be used to select a centroid reference allele.
After this step, the following files and folders should occur in the `analysis`
directory:

* A `GENE` folder, containing modified (merged with flanks) sequences of GENE
  alleles. Sequences are in FASTA format, with `.fa_ex` file extension.
* A `GENE.dist_matrix` file - distance matrix showing the results of
  pairwise alignment for each sequence of the GENE.
* A `GENE.fasta` file - centroid reference sequence of the GENE.
* `database_report` folder, containing the following files:
  * `GENE.lengths` - the distribution of lengths of sequences belonging to
    GENE.
  * `GENE.dist_distr` - the distribution of distances from the aforementioned
    `distance matrix`.
  * `GENE.database_report` - list of similar sequences or subsequences (only
    if similar sequences or subsequences were found).

Once these files are generated, the program will suggest to continue the
calculations. Type "yes" if you want to continue. Type anything else if you
want to stop calculations in order to check the database report files and fix
(delete or add sequences) in the original database. Subsequences or similar
sequences belonging to the same gene can be automatically deleted using command:

    bash clean_database.bh

There might be some additional files (for example, each `GENE.fasta` file is accompanied by
the files, obtained after its indexing).

#### Step 3
If you decided to continue the calculations, then during the third step all
representative sequences will be stored in one `all.fasta` file.
