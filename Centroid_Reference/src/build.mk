
# Makefile for database preprocessing

#SHELL := ./qsub_wrapper

TXT2FA := python ../src/build_database_library.py txt2fa
EXTEND_FASTA := python ../src/build_database_library.py modif_fasta
BUILDMATRIX := python ../src/build_database_library.py dist_matrix
CONSENSUS := python ../src/build_database_library.py compare_consensus
MAKE_REF := python ../src/build_database_library.py make_ref
VCF_MATRIX := python ../src/build_database_library.py vcf_matrix 
COLLECT_OUTLIERS := python ../src/build_database_library.py collect_outliers
CREATEBED := python ../src/build_database_library.py create_bed

INPUT := $(PTH)

#files to remove after all
SAM :=$(addsuffix .sam, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))
BAM :=$(addsuffix .bam, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))
BAM_BAI :=$(addsuffix .sorted.bam.bai, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))
BCF :=$(addsuffix .bcf, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))
VCF_GZ :=$(addsuffix .vcf.gz*, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))

OUTPUT_1 := $(addsuffix .fasta.bwt, $(notdir $(wildcard $(PTH)/*))) $(addsuffix .dist_matrix, $(notdir $(wildcard $(PTH)/*))) $(addsuffix .bed, $(notdir $(wildcard $(PTH)/*)))
OUTPUT_2 := all.fasta.pac $(addsuffix .control, $(notdir $(wildcard $(PTH)/*)))
OUTPUT_3 := $(addsuffix .vcf_matrix, $(notdir $(wildcard $(PTH)/*))) additional_check/outliers

#files to keep
.PRECIOUS: $(addsuffix /dir.done, $(notdir $(wildcard $(PTH)/*)))
.PRECIOUS: $(addsuffix .fa_ex, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa)))) 
.PRECIOUS: $(addsuffix .fasta, $(notdir $(wildcard $(PTH)/*)))
.PRECIOUS: $(addsuffix .fasta.bwt, $(notdir $(wildcard $(PTH)/*)))
.PRECIOUS: $(addsuffix _1.fq, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa)))) 
.PRECIOUS: $(addsuffix _2.fq, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa)))) 
.PRECIOUS: $(addsuffix .vcf, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))
.PRECIOUS: $(addsuffix .consensus, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))
.PRECIOUS: $(addsuffix .sorted.bam, $(basename $(subst $(PTH),, $(wildcard $(PTH)*/*.fa))))

.SECONDEXPANSION:
.SECONDARY:

first_step: $(OUTPUT_1)
second_step: $(OUTPUT_2)
third_step: $(OUTPUT_3)

#removing useless files
clean: clean_bam
	rm $(SAM) 2> /dev/null
clean_bam: clean_bam_bai
	rm $(BAM) 2> /dev/null 
clean_bam_bai: clean_bcf  
	rm $(BAM_BAI) 2> /dev/null
clean_bcf: clean_vcf
	rm $(BCF) 2> /dev/null
clean_vcf: $(PTH)
	rm $(VCF_GZ) 2> /dev/null 
	

#make sure that directories are created
%/dir.done: $(notdir $(wildcard $(PTH)%)) 
	mkdir -p $* && touch $@

#Modify fasta files
%.fa_ex: $(PTH)$$(firstword $$(subst -, , $$*))/$$*.fa $$(firstword $$(subst -, , $$*))/dir.done
	$(EXTEND_FASTA) $< $@ $(FLANKS)

#create distance matrix
%.dist_matrix: $$(addsuffix .fa_ex, $$(basename $$(subst $$(PTH),, $$(wildcard $$(PTH)%/*.fa))))
	$(BUILDMATRIX) "$^" $@ $(PTH) $(FLANKS)

#find the representative sequence
%.fasta: %.dist_matrix 
	$(MAKE_REF) $<

#create the bed-file
%.bed: %.fasta
	$(CREATEBED) $< $@ $(FLANKS)

#index the representative sequence
%.fasta.bwt: %.fasta 
	bwa index $<

#collect all of the representative sequences in one reference sequence
all.fasta: $(addsuffix .fasta, $(notdir $(wildcard $(PTH)*)))
	cat $^ > $@

#index reference sequence
all.fasta.pac: all.fasta
	bwa index $<

#Make fq files
%_1.fq: %.fa_ex
	make_fastq local $(FASTQ) $(basename $<) $< && sleep 10
%_2.fq: %_1.fq
	@
	
#obtain sam file
%.sam: $$(addsuffix .fasta, $$(firstword $$(subst -, , $$*))) %_1.fq %_2.fq
	bwa mem $(BWA) $^ >  $@

#obtain bam file
%.bam:  %.sam
	samtools view -b -S -o $@ $^

#sort bam file
%.sorted.bam: %.bam
	samtools sort $^ -o $@

#index bam file
%.sorted.bam.bai: %.sorted.bam
	samtools index $^

#call variants
%.bcf: $$(addsuffix .fasta, $$(firstword $$(subst -, , $$*))) %.sorted.bam.bai 
	samtools mpileup $(MPILEUP) -g -f $< $(basename $(word 2, $^)) > $@

#save variants in vcf format
%.vcf: %.bcf
	bcftools call -c -v --ploidy $(PLOIDY) $^ > $@

#gzip variants
%.vcf.gz: %.vcf
	bgzip -c $^ > $@

#tabix variants
%.vcf.gz.tbi: %.vcf.gz
	tabix -p vcf $^

#build the consensus seq
%.consensus: $$(addsuffix .fasta, $$(firstword $$(subst -, , $$*))) %.vcf.gz.tbi
	cat $< | vcf-consensus $(basename $(word 2, $^)) > $@

#build the control output
%.control: $$(addsuffix .consensus, $$(basename $$(subst $$(PTH),, $$(wildcard $$(PTH)$$*/*.fa))))
	$(CONSENSUS) "$^" > $@

#build the vcf matrix
%.vcf_matrix: $$(addsuffix .vcf, $$(basename $$(subst $$(PTH),, $$(wildcard $$(PTH)$$*/*.fa)))) %.fasta %.control
	$(VCF_MATRIX) "$^" $@ $(FLANKS)

#collect all outliers into one file
additional_check/outliers: $(addsuffix .control, $(notdir $(wildcard $(PTH)*)))
	mkdir -p additional_check && $(COLLECT_OUTLIERS) "$^" $@





