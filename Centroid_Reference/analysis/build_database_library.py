#!/usr/bin/python

"""
A set of utility scripts to serve the database preprocessing
& analysis.
"""
from __future__ import division
import os
import numpy as np
from fastools import fastools
import argparse
from Bio import Seq, SeqIO, SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

def get_flanks(gene_name, string):
    """
    Returns the flanking regions for the particulat gene
    written in the file. If file containes no information about flanking regions,
    returns poly-N sequence length 50.

    @arg gene_name: the name of gene
    @type gene_name: string
    @arg string: the name of file containing flanking regions
    @type string: string

    """
    all_flanks = open(string).read().split("\n")[:-1]
    flank_l = 'N' * 50
    flank_r = 'N' * 50
    for line in all_flanks:
        fields = line.split()
        if ((fields[0]==gene_name) and (fields[1]=='F')): flank_l = fields[2]
        elif ((fields[0]==gene_name) and (fields[1]=='R')): flank_r = fields[2]
    return flank_l, flank_r
#get_flanks

def modified_fasta(input_handle, output_handle, string):
    """
    Convert FASTA files to the FASTA files
    with flanking regions from both sides. Flanking regions
    should be provided in 'string' file. If no flanking regions provided,
    will be flanked by polyN sequence length 50 from each side.
    Place output FASTA files to the analysis directory,
    returns names of all new files

    @arg input_handle: open readable handle to the input FASTA file
    @type input_handle: stream
    @arg output_handle: open writable handle to the output file
    @type output_handle: stream
    @arg string: name of the file with flanking regions
    @type string: string

    """
    input = SeqIO.read(input_handle, "fasta")
    sequence = input.seq
    id_nmr = input.id
    gene_name = input_handle.name.split(".fasta")[0].split("/")[-1].split("-")[0]
    #gene_name = input_handle.name.split("/")[-2]
    if (string=='NO'):
        flank_l = 'N' * 50
        flank_r = 'N' * 50
    else: flank_l, flank_r = get_flanks(gene_name, string) 
    record = SeqRecord(flank_l+sequence+flank_r, id=id_nmr, description = "") 
    SeqIO.write(record, output_handle, "fasta")
#modified_fasta  

def create_bed(input_handle, output_handle, string):
    """
    Generates bed-file using fasta-file and information about flanking regions.
    Bed-file contain the information about gene, starting and 
    ending point of variable region. Count starts from 0!!

    @arg input_handle: open readable handle to the input FASTA file
    @type input_handle: stream
    @arg output_handle: open writable handle to the output BED file
    @type output_handle: stream
    @arg string: name of the file with flanking regions
    @type string: string
    """
    input = SeqIO.read(input_handle, "fasta")
    sequence_len = len(input.seq)
    gene_name = input_handle.name.split(".fasta")[0].split("/")[-1].split("-")[0]
    if (string=='NO'):
        flank_l = 'N' * 50
        flank_r = 'N' * 50
    else: flank_l, flank_r = get_flanks(gene_name, string)
    coord_l = len(flank_l)
    coord_r = sequence_len - 1 - len(flank_r)
    output_handle.write("{0}\t{1}\t{2}".format(gene_name, coord_l, coord_r)) 
#create_bed

def dist_matrix(input_files, output_handle, string, string_2):
    """
    Build the matrix with paiwise Levenshtein distance
    between all FASTA files in the provided list. 
    Saves result to the output file. Also creates
    the report about input sequences (quality, length
    and distance distributions)
    
    @arg input_files: list of file names
    @type input_files: str
    @arg output_handle: open writable handle to the output file
    @type output_handle: stream
    @arg string: path to the original sequences
    @type string: string
    @arg string_2: name of the file with flanking regions
    @type string_2: string

    """
    for_check= []
    dir = input_files.split()[0].split("/")[-2]
    lengths = open("database_report/{0}.lengths".format(dir), "w")
    distances = []
    if (string_2=='NO'):
        flank_l = 'N' * 50
        flank_r = 'N' * 50
    else: flank_l, flank_r = get_flanks(dir, string_2)

    for name_1 in input_files.split():
        output_handle.write("{0}".format(name_1.split(".fa_ex")[0]))
        seq1 = SeqIO.read(name_1, "fasta").seq[len(flank_l):-(len(flank_r))]

        for name_2 in input_files.split():
            if (name_1 == name_2): output_handle.write("\t0")
            else:
                seq2 = SeqIO.read(name_2, "fasta").seq[len(flank_l):-(len(flank_r))]
                check = [name_1, name_2]
                rev_check = [name_2, name_1]
                distance = fastools.aln(check)[0][2]
                distances.append(distance)
                if ((check not in for_check) and 
                    ((seq1 in seq2) or (seq2 in seq1) or distance == 0)):
                        report_subsequence(check, [seq1, seq2], dir, string)
                        for_check.append(check)
                        for_check.append(rev_check)
                output_handle.write("\t{0}".format(distance))

        output_handle.write("\n")
        lengths.write("{0}\t{1}\n".format(name_1, len(seq1)))
    
    if (len(input_files.split())>=5):
        distrib = open("database_report/{0}.dist_distr".format(dir), "w")
        hist = np.histogram(distances, bins = 20)
        for position in range(len(hist[0])): 
            distrib.write("{0}\t{1}\n".format(hist[1][position], hist[0][position]/2))
        distrib.close()
    
    lengths.close()
#dist_matrix

def vcf_matrix(input_files, output_handle, string):
    """
    Build the matrix with paiwise  distance
    between all vcf files  (except the ones in the file with outliers) 
    in the provided list.
    Saves result to the output file.

    @arg input_files: list of file names
    @type input_files: str
    @arg output_handle: open writable handle to the output file
    @type path: stream
    @arg string: name of the file with flanking regions
    @type string: string

    """
    outliers = open(input_files.split()[-1]).read().split("\n")[:-2]
    len_total = len(SeqIO.read(input_files.split()[-2], "fasta").seq)
    gene_name = input_files.split()[-2].split(".fasta")[0].split("/")[-1].split("-")[0]
    print (gene_name)
    if (string=='NO'):
        flank_l = 'N' * 50
        flank_r = 'N' * 50
    else: flank_l, flank_r = get_flanks(gene_name, string)
    len_l = len(flank_l)
    len_r = len(flank_r)
    for name_1 in input_files.split()[:-2]:
        if (name_1.split(".")[0] not in outliers):
            line_header = name_1.split(".")[0]
            output_handle.write("{0}".format(line_header))
            for name_2 in input_files.split()[:-2]:
                if (name_2.split(".")[0] not in outliers):
                    distance = vcf_distance_self("{0} {1} {2} {3} {4}".format(name_1, name_2, len_l, len_r, len_total))
                    output_handle.write("\t{0}".format(distance))
            output_handle.write("\n")
#vcf_matrix

def analyse_dist_matrix(input_handle):
    """
    Reports the allele with the lowest average alignment score
    witheen all sequences in matrix, record the corresponding 
    sequence with gene ID to the fasta file. Check the distance 
    distributions within the matrix file as well 

    @arg input_handel: open readable handle to a matrix file
    @type files: stream

    """
    result = {}
    lines = input_handle.read().split("\n")[:-1]
    name = input_handle.name.split(".")[0]
    for line in lines:
        scores = []
        all = line.split("\t")
        for i in range(1,len(all)):
            scores.append(int(line.split("\t")[i]))
        av = np.mean(np.array(scores))
        result[line.split("\t")[0]] = av
    file = sorted(result.items(), key=lambda result: result[1])[0][0] + ".fa_ex"
    seq = SeqIO.read(file, "fasta").seq
    record = SeqRecord(seq, name, "", "")
    SeqIO.write(record, name + ".fasta", "fasta")
# analyse_dist_matrix

def analyse_vcf(input_files, string):
    """
    Compare the given vcf file with every vcf file 
    from the given matrix. Print out the closest vcf file +
    the average distance to the closest neighbour within the 
    matrix file

    @arg input_files: vcf file and matrix file
    @type inpit_files: str
    @arg string: name of the file with flanking regions
    @type string: str

    """
    matrix, file = input_files.split()
    average = av_distance_closest(matrix) 

    gene_name = matrix.split(".")[0]
    len_total = len(SeqIO.read(gene_name + ".fasta", "fasta").seq)
    if (string=='NO'):
        flank_l = 'N' * 50
        flank_r = 'N' * 50
    else: flank_l, flank_r = get_flanks(gene_name, string)
    len_l = len(flank_l)
    len_r = len(flank_r)

    closest, distance = closest_vcf(matrix, file, len_l, len_r, len_total)
    heterozygous = check_zigosity(file, len_l, len_r, len_total) 
    if (heterozygous==0): print ("| {0} | {1} | {2} |".format(closest, distance, average))
    if (heterozygous>0): print ("| {0} | {1} | {2} | {3} NON-HOMOZYGOUS SITES DETECTED!!! |".format(closest, distance, average, heterozygous))
#analyse_vcf

def check_zigosity(file, len_l, len_r, len_total):
    """
    Checks whether or not there are any signs of heterozygosity in the 
    vcf-file

    @arg file: vcf file
    @type file: string
    @arg len_l: lengh of the left flank
    @type len_l: int
    @arg len_r: lengh of the right flank
    @type len_r: int
    @arg len_total: lengh of the reference 
    @type len_total: int
    @type returns: int

    """
    heterozygosity = 0
    lines = open(file).read().split(".bam")[-1].split("\n")[1:-1]
    for line in lines:
        line_split = line.split('\t')
        if (int(line_split[1])>int(len_l) and
            int(line_split[1])<=int(len_total)-int(len_r) and 
            ("AF1=1" not in line)): heterozygosity+=1
    return heterozygosity
#check_zygosity

def coverage(input_files, string, string_2):
    """
    Returns the horisontal coverage using the wig file:
    ratio between the number of positions with coverage
    >= string value and total number of positions

    @arg input_files: wig file and reference file
    @type files: str
    @arg string: threshold for the coverage
    @type string: str
    @arg string_2: the files with flanks
    @type string_2: str

    """
    positions = {}
    wig, reference = input_files.split()
    seq = SeqIO.read(reference, "fasta").seq
    length = len(seq)
    
    gene_name = reference.split(".fasta")[0].split("/")[-1].split("-")[0]
    #gene_name = reference.split(".")[0]
    if (string_2=='NO'):
        flank_l = 'N' * 50
        flank_r = 'N' * 50
    else: flank_l, flank_r = get_flanks(gene_name, string_2)
    len_l = len(flank_l)
    len_r = len(flank_r)
    data = open(wig).read().split("\n")[2:-1]
    for i in data: positions[i.split()[0]]=i.split()[1]
    k = 0
    for position in range(len_l+1, length-len_r+1):
        if (str(position) in positions.keys()):
            if (int(positions[str(position)]) >= int(string)): k+=1
    print (round(k/(length - len_l - len_r), 3))
# coverage

def list_for_removal(input_handle):
    """
    Processing the database repport file and suggest
    files to be removed from the database

    @arg input_handle: open readable handle to the database repport file
    @type input_handle: stream

    """
    for_del = []
    lines = input_handle.read().split("\n")[:-1]
    name = input_handle.name.split("_")[-1]
    for i in lines:
        if ("are equal" in i): rem = i.split(" and ")[-1].split()[0]
        elif("is a subsequence" in i): rem = i.split(" is a subsequence")[0].split()[-1]
        if (rem not in for_del): for_del.append(rem)

    out = '{0}.remove'.format(name)
    f = open(out, "w")
    for delet in for_del: f.write("{0}\n".format(delet))
    f.close()
    print (out)
#list_for_removal

def report_subsequence(list_files, list_sequences, dir_name, path):
    """
    reports whether 2 provided fasta files have the same content
    or whether one of the sequences is a subsequence of another one
    
    @arg list_files: list with 2 fasta files
    @type list_files: list
    @arg list_sequences: list with 2 fasta files
    @type list_sequences: list with corresponding sequences
    @arg dir_name: name of the gene, for which the report file
                   will be created
    @type dir_name: string
    @arg path: path to the original database
    @type path: string

    """
    f = open("database_report/{0}.database_report".format(dir_name), "a")
    if (len(list_sequences[0]) == len(list_sequences[1])):
        f.write("Sequences {2}{0} and {2}{1} are equal\n".format(list_files[0].split("_ex")[0], list_files[1].split("_ex")[0], path))
    elif (len(list_sequences[0]) > len(list_sequences[1])):
        f.write("Sequence {2}{1} is a subsequence of {2}{0}\n".format(list_files[0].split("_ex")[0], list_files[1].split("_ex")[0], path))
    elif (len(list_sequences[0]) < len(list_sequences[1])):
        f.write("Sequence {2}{0} is a subsequence of {2}{1}\n".format(list_files[0].split("_ex")[0], list_files[1].split("_ex")[0], path))
    f.close()
#report_subsequence

def av_distance(matrix):
    """
    calculates the average distance in a square
    distance matrix

    @arg matrix: name of the input file
    @type matrix: string
    @returns: average distance in matrix
    @type returns: float

    """
    
    sum = 0
    lines = open(matrix).read().split("\n")[:-1]
    for line in lines:
        for value in line.split()[1:]: sum = sum + float(value)
    return round((sum * 2)/(len(lines) * (len(lines) - 1)), 5)
#av_distance

def av_distance_closest(matrix):
    """
    calculates the average distance to the closest neighbour
    in a square distance matrix

    @arg matrix: name of the input file
    @type matrix: string
    @returns: average distance to the closest neighbour in matrix
    @type returns: float

    """

    sum = 0
    lines = open(matrix).read().split("\n")[:-1]
    for line in lines:
            all = []
            for value in line.split()[1:]: all.append(float(value))
            sum = sum + sorted(all)[1]
    return round(sum/len(lines), 5)
#av_distance_closest

def vcf_distance_self(input_files):
    """
    Calculates ammounts of mismatches between 2 given vcf files

    @arg input_files: list with vcf files names, lengths of flanks and total length
    @type input_files: list
    @type returns: int

    """
    list_ = input_files.split()
    mismatch = 0
    if (list_[0] == list_[1]):
        mismatch = 0
    else:
        counts_file_1 = count_vcf_variants(list_[0], list_[2], list_[3], list_[4])
        counts_file_2 = count_vcf_variants(list_[1], list_[2], list_[3], list_[4])
        keys_1 = set(counts_file_1.keys())
        keys_2 = set(counts_file_2.keys())
        intersection = keys_1 & keys_2
        if (len(intersection)==0): mismatch = len(keys_1) + len(keys_2)
        else:
            similar = 0
            for value in intersection:
                if (counts_file_1[value]==counts_file_2[value]): similar+=2
            mismatch = len(keys_1) + len(keys_2) - similar
    return(mismatch)
#vcf_distance_self

def compare_consensus(input_files):
    """
    Compares each consensus file from the list with the original fasta file
    Prints out names for which sequences doesn't match

    @arg vcf_list: list with consensus files names
    @type vcf_list: list

    """
    list_ = input_files.split()
    for line in list_:
        name = line.split(".")[0]
        fa_name = ".".join([name, "fa_ex"])
        record1 = SeqIO.read(line, "fasta").seq
        record2 = SeqIO.read(fa_name, "fasta").seq
        if (str(record1)!=str(record2)): print (name)
    print ("Analysis is ended")
#compare_consensus

def count_vcf_variants(file_name, len_l, len_r, len_total):
    """
    Counts amount of varints in the provided vcf file.
    Returns the dictionary with variant and position

    @arg file_name: name of the input vcf-file
    @type file_name: string
    @arg len_l: length of the left flank
    @type len_l: int
    @arg len_r: length of the right flank
    @type len_r: int
    @arg len_total: length of the reference
    @type len_total: int

    """
    count = {}
    empty_or_not = open(file_name).read().split('.bam')
    if (len(empty_or_not)>1):
        lines = empty_or_not[-1].split('\n')[1:-1]
        for line in lines:
            line_split = line.split('\t')
            if (int(line_split[1])>int(len_l) and 
                int(line_split[1])<=int(len_total)-int(len_r)): count[int(line_split[1])] = (line_split[4].split(',')[0])
    return count
#count_vcf_variants

def count_vcf(input_files, string):
    """
    Counts amount of varints in the provided vcf file.
    Ignores 'X' symbol in ALT column. Returns the ammounts
    of variants in the vcf-file

    @arg input_files: vcf file and reference file
    @type inpit_files: str
    @arg string: name of the file with flanking regions
    @type string: str

    """
    file, reference = input_files.split()
    gene_name = reference.split(".fasta")[0].split("/")[-1].split("-")[0]
    flanks = string
    if (flanks=='NO'):
        flank_l = 'N' * 50
        flank_r = 'N' * 50
    else: flank_l, flank_r = get_flanks(gene_name, flanks)
    len_l = len(flank_l)
    len_r = len(flank_r)
    len_total = len(SeqIO.read(reference, "fasta").seq)
    print (len(count_vcf_variants(file, len_l, len_r, len_total)))
#count_vcf

    
def closest_vcf(matrix, file, len_l, len_r, len_total):
    """
    Compares the given vcf file with every vcf file from the matrx,
    returns the name of the closest file

    @arg matrix: matrix file
    @type matrix: string
    @arg file: vcf file
    @type file: string
    @arg len_l: lengh of the left flank
    @type len_l: int
    @arg len_r: lengh of the right flank
    @type len_r: int
    @arg len_total: lengh of the reference
    @type len_total: int
    @type returns: variables

    """
    names = []
    values = []
    lines = open(matrix).read().split("\n")[:-1]
    for line in lines:
        name = line.split("\t")[0]
        names.append(name)
        values.append(vcf_distance_self("{0} {1} {2} {3} {4}".format(name+".vcf", file, len_l, len_r, len_total)))
    index = values.index(min(values))
    return names[index], values[index]
#closest_vcf

def collect_outliers(input_files, output_handle):
    """
    Collects the informationa about outliers from the 
    input files. Sequences of the outliers are stored on the 
    "additional_check" dirrectory together with output files
    which is the list of all outliers

    @arg input_files: list of control files
    @type input_files: string
    @arg output_handle: open writable handle to the output file
    @type path: stream

    """
    for name in input_files.split():
        list_of_files = open(name).read().split("\n")[:-1]
        for i in list_of_files[:-1]: 
            output_handle.write("{0}\n".format(i))
            output_name = "additional_check/"+i.split("/")[-1]+".fasta"
            record = SeqIO.read(i+".fa_ex", "fasta")
            SeqIO.write(record, output_name, "fasta")            
#vcf_matrix

def compare_analysis_parts(input_files, output_handle):
    """
    Compare results of analysis for the vcf-independent and vcf-dependent 
    parts of the analysis. Saves results on the comparison to the outpput file

    @arg input_files: list of files (outcom of vcf-dependent and vcf-independent parts)
    @type input_files: string
    @arg output_handle: open writable handle to the output file
    @type path: stream

    """
    main_file, additional_file = input_files.split()
    main = open(main_file).read().split("\n")[:-1]
    
    main_parsed = {}
    for k in main:
        values_ma = k.split(" | ")
        if ("IS NOT FOUND" in k): main_parsed[values_ma[1]]="IS NOT FOUND" 
        else: main_parsed[values_ma[0]]=[float(values_ma[1]), values_ma[2], int(values_ma[3]), float(values_ma[4].split(" |")[0])]

    additional = open(additional_file).read().split("\n")[:-1]
    for i in additional:
        values_ad = i.split()
        name = values_ad[0]
        gene = name.split("/")[-1].split("-")[0]
        cov = float(values_ad[1])
        mism = int(values_ad[2])
        if (gene not in main_parsed.keys()): main_parsed[gene]=[cov, name, mism, "from the low-homology group"]
        elif (main_parsed[gene]=="IS NOT FOUND"): main_parsed[gene]=[cov, name, mism, "from the low-homology group"]
        else:
            if (cov>=main_parsed[gene][0] and mism<main_parsed[gene][2]): main_parsed[gene]=[cov, name, mism, "from the low-homology group"]

    for g in main_parsed: output_handle.write("{0} | {1} | {2} | {3} | {4} |\n".format(g, main_parsed[g][0], main_parsed[g][1], main_parsed[g][2], main_parsed[g][3]))
#compare_analysis_parts

def main():
    """
    Command line interface.
    :arg args: Arguments passed to :meth:`argparse.ArgumentParser.parse_args`
    (default is to use `sys.argv[1:]`).
    :type args: list(str)

    """

    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument('input_handle', metavar='INPUT', type=argparse.FileType('r'),
                              help='input file')

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument('output_handle', metavar='OUTPUT', type=argparse.FileType('w'),
                               help='output file')

    input_files_parser = argparse.ArgumentParser(add_help=False)
    input_files_parser.add_argument('input_files', metavar='LIST', type=str,
                                    help='list of input files')

    str_parser = argparse.ArgumentParser(add_help=False)
    str_parser.add_argument('string', metavar='STR', type=str,
                            help='string')
    str_parser_2 = argparse.ArgumentParser(add_help=False)
    str_parser_2.add_argument('string_2', metavar='STR', type=str,
                              help='string_2')

    list_parser = argparse.ArgumentParser(add_help=False)
    list_parser.add_argument('list', metavar='LIST', type=list,
                             help='list')

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers()

    parser_analyse_dist_matrix = subparsers.add_parser('make_ref',
                                                       parents=[input_parser])
    parser_analyse_dist_matrix.set_defaults(func=analyse_dist_matrix)

    parser_analyse_vcf = subparsers.add_parser('analyse_vcf',
                                               parents=[input_files_parser, str_parser])
    parser_analyse_vcf.set_defaults(func=analyse_vcf)

    parser_collect_outliers = subparsers.add_parser('collect_outliers',
                                                    parents=[input_files_parser, output_parser])
    parser_collect_outliers.set_defaults(func=collect_outliers)

    parser_compare_analysis_parts = subparsers.add_parser('compare_analysis_parts',
                                                          parents=[input_files_parser, output_parser])
    parser_compare_analysis_parts.set_defaults(func=compare_analysis_parts)

    parser_compare_consensus = subparsers.add_parser('compare_consensus',
                                                     parents=[input_files_parser])
    parser_compare_consensus.set_defaults(func=compare_consensus)

    parser_coverage = subparsers.add_parser('coverage',
                                            parents=[input_files_parser, str_parser, str_parser_2])
    parser_coverage.set_defaults(func=coverage)
    
    parser_count_vcf = subparsers.add_parser('count_vcf',
                                             parents=[input_files_parser, str_parser])
    parser_count_vcf.set_defaults(func=count_vcf)

    parser_create_bed = subparsers.add_parser('create_bed',
                                              parents=[input_parser, output_parser, str_parser])
    parser_create_bed.set_defaults(func=create_bed)
    
    parser_dist_matrix = subparsers.add_parser('dist_matrix',
                                               parents=[input_files_parser, output_parser, str_parser, str_parser_2])
    parser_dist_matrix.set_defaults(func=dist_matrix)

    parser_list_for_removal = subparsers.add_parser('rem_list',
                                                     parents=[input_parser])
    parser_list_for_removal.set_defaults(func=list_for_removal)

    parser_modified_fasta = subparsers.add_parser('modif_fasta',
                                                  parents=[input_parser, output_parser, str_parser])
    parser_modified_fasta.set_defaults(func=modified_fasta)

    parser_vcf_dist_self = subparsers.add_parser('vcf_dist_self',
                                                 parents=[input_files_parser])
    parser_vcf_dist_self.set_defaults(func=vcf_distance_self)

    parser_vcf_matrix = subparsers.add_parser('vcf_matrix',
                                              parents=[input_files_parser, output_parser, str_parser])
    parser_vcf_matrix.set_defaults(func=vcf_matrix)


    try:
        arguments = parser.parse_args()
    except IOError as error:
        parser.error(error)

    try:
        arguments.func(**dict((k, v) for k, v in vars(arguments).items()
                              if k not in ('func', 'subcommand')))
    except ValueError as error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
