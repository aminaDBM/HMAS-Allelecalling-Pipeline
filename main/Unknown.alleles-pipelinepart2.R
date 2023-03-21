
#!/usr/bin/env Rscript

library(conflicted)
library(vegan)
library(tidyverse)
library(dplyr)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::lag)
library(stringr)
library(openxlsx)
library("seqinr")

#setwd("Z:\\266Pacbio.IlluminaData\\data\\unknown_alleles\\analysis\\")


#creating lists of input files 
tsv_files <- list.files(pattern='*.unknown_Alleles.tsv')
uniqseq_files <- list.files(pattern='*.unknown_Alleles.fasta')


#initiating list variables
ListofNovelAlleles <- list()
no.ofUnknownAlleles <- list()
no.ofNovelAlleles <- list()
no.ofAllelesNotfound <- list()



for(i in 1:length(tsv_files)) {
  prefix <- str_remove(tsv_files, pattern = ".unknown_Alleles.tsv")
  
  ################## READ IN DATA ################# 
  
  #BLAST result from CentroidRefDB
  
  blastoutput <- read.delim(file = tsv_files[i],
                             header = FALSE,
                             sep = "\t",
                             stringsAsFactors = FALSE)
  
  
  #Read in fasta file *.unknownAlleles.fasta
  seqfile <- read.fasta(file = uniqseq_files[i], forceDNAtolower = FALSE) 
  
  
  ##################################################
  #saving ".tsv" files as a dataframe "blastoutput"
  
  #naming the cols
  colnames(blastoutput) <- c('qseqid', 'sseqid', 'pident', 'qcovs', 'no.ofgaps', 'BTOP', 'qseqlen', 'sseqlen', 'alignlen')
  
  #calculating %age length
  blastoutput$percentage_length <- ((blastoutput$alignlen - blastoutput$no.ofgaps)/(blastoutput$sseqlen-100))*100 #making a col for sample name
  
  #Applying filters to identify novel alleles
  novel_alleles <- dplyr::filter(blastoutput, str_detect(qseqid, sseqid)) #checking the primerpair match
  excluded_alleles <- dplyr::filter(blastoutput, !( str_detect(qseqid, sseqid))) #checking the primerpair match
  
  novel_alleles <- novel_alleles %>% filter(pident >= 85) #checking homology >=85
  novel_alleles <- novel_alleles %>% filter(percentage_length >= 90) #checking %age length >=90
  novel_alleles <- novel_alleles %>% filter(no.ofgaps <= 10) #checking no.of gaps <=10
  
  #writing dataframe to csv file and appending with each run
  write.csv(novel_alleles, file = paste0(prefix[i],".novelAlleleStatistics.csv"), row.names=FALSE)

  ############################################################
  
  
  #Annontating fasta files and writing output fasta
  
  #identifying the list of novel alleles
  ListofNovelAlleles[[i]] <- unique(novel_alleles$qseqid)
  seqID <- unlist(ListofNovelAlleles[i])

  
  #Splitting the unknown alleles into novel and not_found alleles
  
  novel_Alleles <- seqfile[names(seqfile) %in% seqID]
  notfound_Alleles <- seqfile[!names(seqfile) %in% seqID]

  #Calculating No. of sequences in : UnknownAlleles(INPUT), NovelAlleles(OUTPUT1), AllelesNotFound(OUTPUT2)
  no.ofUnknownAlleles[[i]] <- length(seqfile)
  no.ofNovelAlleles[[i]] <- length(novel_Alleles)
  no.ofAllelesNotfound[[i]] <- length(notfound_Alleles)
  
  #writing fasta file for novel_alleles.fasta and alleles NotFound
  write.fasta(sequences = novel_Alleles,names=names(novel_Alleles),file.out = paste0(prefix[i],".novel_Alleles.fasta"))
  
  if (length(notfound_Alleles) != 0){
    write.fasta(sequences=notfound_Alleles,names=names(notfound_Alleles),file.out = paste0(prefix[i],".notfound_Alleles.fasta"))
    }
  
  
  #writing output from all input samples into one consolidated output file
  nestedlist <- list(prefix, no.ofUnknownAlleles, no.ofNovelAlleles, no.ofAllelesNotfound)
  OutputStats <- as.data.frame(do.call(cbind, nestedlist))
  #naming the cols
  colnames(OutputStats) <- c('Sample', 'Unknown Alleles', 'Novel Alleles', 'Not Found')
  OutputStats
  
}                                    

df_to_csv <- data.frame(lapply(OutputStats, as.character), stringsAsFactors=FALSE)  
write.csv(df_to_csv, file = "UnknownAlleles.consolidatedStats.csv", row.names=FALSE)
