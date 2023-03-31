
#!/usr/bin/env Rscript
#args=(commandArgs(TRUE))

#path<-Sys.getenv("Path_to_Input")


library(conflicted)
library('seqinr')
library(vegan)
library(tidyverse)
library(dplyr)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::lag)

setwd(getwd())

#creating lists of input files 
tsv_files <- list.files(pattern='_extractedAmplicons.tsv')
uniqseq_files <- list.files(pattern='_extractedAmplicons.fasta')


#initiating list variables
ListofKnownAlleles <- list()
no.ofunknownAlleles <- list()
no.ofknownAlleles <- list()
no.ofsampleSeqs <- list()

for(i in 1:length(tsv_files)) {                                     # Head of for-loop
  
  prefix <- str_remove(tsv_files, pattern = ".tsv")
  print(tsv_files[i])
  #BLAST result with pident=100
  blastoutput1 <- read.delim(file = tsv_files[i],
                             header = FALSE,
                             sep = "\t",
                             stringsAsFactors = FALSE)
  
  ####################################################
  
  #naming the cols and arranging table by qcovs values
  colnames(blastoutput1) <- c('qseqid', 'sseqid', 'pident', 'qcovs', 'qseqlen', 'sseqlen', 'alignlen')
  blastoutput1[c('PrimerGroup', 'Isolate')] <- str_split_fixed(blastoutput1$sseqid, '-', 2)#making a col for sample name
  
  #Filtering known alleles: on basis of qseqlen=sseqlen=alignlen to get the sequences which align with the unique insilico primer
  #sequence 
  blastoutput1$All_equal <- blastoutput1$qseqlen==blastoutput1$sseqlen & blastoutput1$sseqlen==blastoutput1$alignlen
  known_alleles <- filter (blastoutput1, blastoutput1$All_equal=='TRUE') #qseqlen=sseqlen=alignlen
  known_alleles <- dplyr::filter(known_alleles, grepl('100', qcovs)) #keeping only the qcovs=100 values
  known_alleles <- dplyr::filter(known_alleles, str_detect(qseqid, PrimerGroup)) #checking the primerpair match
  ListofKnownAlleles[[i]] <- unique(known_alleles$qseqid)
  
  #writing dataframe to csv file and appending with each run
  write.csv(known_alleles, file = paste0(prefix[i],"-KnownAlleleStatistics.csv"),  row.names = FALSE)
  ListofKnownAlleles
#######################################################

  seqID <- unlist(ListofKnownAlleles[i])
  #Unique.seqs (mothur high quality data file)
  seqfile = read.fasta(file = uniqseq_files[i], forceDNAtolower = FALSE)
  
  #Splitting the unique.fasta into known and unknown alleles
  Known_fasta <- seqfile[names(seqfile) %in% seqID]
  Unknown_fasta <- seqfile[!names(seqfile) %in% seqID]

  #Calculating No. of: UnknownAlleles(INPUT), NovelAlleles(OUTPUT1), AllelesNotFound(OUTPUT2)
  no.ofsampleSeqs[[i]] <- length(seqfile)
  no.ofknownAlleles[[i]] <- length(Known_fasta)
  no.ofunknownAlleles[[i]] <- length(Unknown_fasta)
  
  #writing fasta file for known_alleles.fasta and unknown_alleles.fasta
  write.fasta(sequences=Known_fasta,names=names(Known_fasta),file.out=paste0(prefix[i],".known_Alleles.fasta"))

  #writing unknown_allele.fasta only when the list is not NULL
  if(no.ofunknownAlleles[i] > 0){
    
    write.fasta(sequences=Unknown_fasta,names=names(Unknown_fasta),file.out=paste0(prefix[i],".unknown_Alleles.fasta"))
  }
  
  
  
  #writing a consolidated csv file for all the samples
  nestedlist <- list(prefix, no.ofsampleSeqs, no.ofknownAlleles, no.ofunknownAlleles)
  OutputStats <- as.data.frame(do.call(cbind, nestedlist))
  #naming the cols
  colnames(OutputStats) <- c('Sample', 'No.of Sample sequences', 'Known Alleles', 'Unknown Alleles')
  OutputStats
}
  
df_to_csv <- data.frame(lapply(OutputStats, as.character), stringsAsFactors=FALSE)  
write.csv(df_to_csv, file = "KnownAlleles.consolidatedStats.csv", row.names=FALSE)  

   
  

  








