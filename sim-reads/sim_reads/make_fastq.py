#!/usr/local/bin/python

# TODO Make the quality score a command line parameter.

import argparse
import random
import suds

from Bio import Seq, SeqIO

from . import doc_split, usage, version


# Location of the Mutalyzer webservice.
mutalyzer_service_description = 'https://mutalyzer.nl/services/?wsdl'
read_number = 1


def _c2r(reference_length, coverage, read_length):
    """
    Calculate the number of fragments based on the desired coverage.

    :arg int reference_length: Lenght of the reference sequence.
    :arg int coverage: Desired coverage.
    :arg int read_length: Size of the reads.

    :return int: Number of fragments.
    """
    return (reference_length * coverage) // (read_length * 2)


def _get_variant(handle):
    """
    Read one line from the input file, strip the newline and return it.

    :arg stream handle: Open handle to the input file.
    """
    return handle.readline().strip('\n')


def _write_fastq(results, reference, number_of_fragments, insert_size,
        variance, read_length):
    """
    Make simulated reads.

    :arg list[stream] results: List of open writable file handles.
    :arg str reference: The reference sequence.
    :arg int number_of_fragments: Number of read pairs to simulate.
    :arg int insert_size: Approximate distance between the reads.
    :arg int variance: Variation of the read length.
    :arg int read_length: Size of the reads.
    """
    global read_number

    for i in range(number_of_fragments):
        this_insert_size = max(read_length,
            int(random.normalvariate(insert_size, variance)))
        position = random.randint(0, len(reference) - this_insert_size)

        results[0].write('@%i/1\n%s\n+\n%s\n' % (read_number, 
            reference[position:position + read_length], 'b' * read_length))
        results[1].write('@%i/2\n%s\n+\n%s\n' % (
            read_number, Seq.reverse_complement(str(
                reference[position + this_insert_size - read_length:
                    position + this_insert_size])), 'b' * read_length))
        read_number += 1


def mutate(input_handle, output, insert_size, variance, read_length,
        number_of_fragments, coverage, reference, start, end, accno,
        orientation, heterozygous):
    """
    Use a file containing variants to obtain a mutated reference sequence from
    Mutalyzer. Then make paired end reads out of the mutated sequence.

    :arg stream input_handle: Open readable handle to the file containing
        variant descriptions.
    :arg str output: Prefix of the names of the output files.
    :arg int insert_size: Approximate distance between the reads.
    :arg int variance: Variation of the read length.
    :arg int read_length: Size of the reads.
    :arg int number_of_fragments: Number of read pairs to simulate.
    :arg int coverage: Desired coverage.
    :arg str reference: Chromosomal accession number.
    :arg int start: Begin of the slice of the reference sequence.
    :arg int end: End of the slice of the reference sequence.
    :arg str accno: Custom accession number.
    :arg int orientation: Orientation of the slice (1=forward, 2=reverse).
    :arg bool heterozygous: Make heterozygous variants.
    """
    # Read the input file.
    allelic_variant = _get_variant(input_handle)
    line = _get_variant(input_handle)
    while line:
        allelic_variant = '%s;%s' % (allelic_variant, line)
        line = _get_variant(input_handle)
    
    # Set up the SOAP interface to Mutalyzer.
    mutalyzer_service = suds.client.Client(
        mutalyzer_service_description).service

    # Retrieve the reference sequence.
    if not accno:
        accno = mutalyzer_service.sliceChromosome(chromAccNo=reference,
            start=start, end=end, orientation=orientation)

    # Mutate the reference sequence.
    mutalyzer_output = mutalyzer_service.runMutalyzer(
        variant = '%s:g.[%s]' % (accno, allelic_variant))
    sequence = mutalyzer_output.mutated

    # Write the chromosomal description to the results file.
    results_handle = open('%s.txt' % output, 'w')
    results_handle.write('%s: %s:%i_%i\n\n' % (accno, reference, start, end))
    
    if int(mutalyzer_output.errors) > 0:
        for i in mutalyzer_output.messages[0]:
            raise ValueError('(%s): %s' % (i.errorcode, i.message))
    if 'chromDescription' in mutalyzer_output:
        description = mutalyzer_output.chromDescription
    else:
        description = mutalyzer_output.genomicDescription

    for i in description.split('[')[1][:-1].split(';'):
        results_handle.write('%s\n' % i)
    results_handle.close()

    results = [open('%s_1.fq' % output, 'w'), open('%s_2.fq' % output, 'w')]
    number_of_pairs = (_c2r(len(sequence), coverage, read_length) or
        number_of_fragments)
    if heterozygous:
        _write_fastq(results, sequence, number_of_pairs / 2, insert_size,
            variance, read_length)
        _write_fastq(results, mutalyzer_output.original, number_of_pairs / 2,
            insert_size, variance, read_length)
    else:
        _write_fastq(results, sequence, number_of_pairs, insert_size, variance,
            read_length)


def local(output, reference_handle, insert_size, variance, read_length,
        number_of_fragments, coverage):
    """
    Use a local fasta file to make paired end reads.

    :arg str output: Prefix of the names of the output files.
    :arg stream reference_handle: Open readable handle to the reference file.
    :arg int insert_size: Approximate distance between the reads.
    :arg int variance: Variation of the read length.
    :arg int read_length: Size of the reads.
    :arg int number_of_fragments: Number of read pairs to simulate.
    :arg int coverage: Desired coverage.
    """
    results = [open('%s_1.fq' % output, 'w'), open('%s_2.fq' % output, 'w')]

    for record in SeqIO.parse(reference_handle, 'fasta'):
        number_of_pairs = (_c2r(len(record.seq), coverage, read_length) or
            number_of_fragments)
        _write_fastq(results, str(record.seq), number_of_pairs, insert_size,
            variance, read_length)


def main():
    """
    Main entry point.
    """
    parent_parser = argparse.ArgumentParser('parent', add_help=False)
    parent_parser.add_argument('output', metavar='OUTPUT', type=str,
        help='prefix of the names of the output files')
    parent_parser.add_argument('-s', dest='insert_size', type=int, default=300,
        help='mean insert size (default=%(default)s)')
    parent_parser.add_argument('-d', dest='variance', type=int, default=25,
        help='standard deviation of the insert size (default=%(default)s)')
    parent_parser.add_argument('-l', dest='read_length', type=int, default=50,
        help='read length (default=%(default)s)')
    parent_parser.add_argument('-n', dest='number_of_fragments', type=int,
        default=1000000, help='number of fragments (default=%(default)s)')
    parent_parser.add_argument('-c', dest='coverage', type=int,
        default=0, help='coverage (default=%(default)s)')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action='version', version=version(parser.prog))
    subparsers = parser.add_subparsers()

    parser_mutate = subparsers.add_parser('mutate', parents=[parent_parser],
        description=doc_split(mutate), epilog="""If a chromosomal accession
        number (option -r) is used, the options -b and -e are used to make a
        slice of this chromosome. To retrieve the reference sequence of a gene,
        use the -u option (the options -r, -b and -e are ignored in this case).
        With the OUTPUT option the prefix for three output files is given:
        OUTPUT.txt, OUTPUT_1.fq and OUTPUT_2.fq""")
    parser_mutate.add_argument('input_handle', metavar='INPUT',
        type=argparse.FileType('r'), help='name of the input file')
    parser_mutate.add_argument('-r', dest='reference', default='NC_000008.10',
        type=str, help='chromosomal accession number (default=%(default)s)')
    parser_mutate.add_argument('-b', dest='start', type=int, default=136800000,
         help='begin of the slice (default=%(default)s)')
    parser_mutate.add_argument('-e', dest='end', type=int, default=139000000,
        help='end of the slice (default=%(default)s)')
    parser_mutate.add_argument('-u', dest='accno', type=str, default='',
         help='accession number of a referencence sequence')
    parser_mutate.add_argument('-t', dest='orientation', default=1, const=2,
        action='store_const', help='reverse the orientation of the slice')
    parser_mutate.add_argument('-H', dest='heterozygous', default=False,
        action='store_true', help='make heterozygous variants')
    parser_mutate.set_defaults(func=mutate)

    parser_local = subparsers.add_parser('local', parents=[parent_parser],
        description=doc_split(local))
    parser_local.add_argument('reference_handle', metavar='REFERENCE',
        type=argparse.FileType('r'), help='name of a local reference sequence')
    parser_local.set_defaults(func=local)

    arguments = parser.parse_args()

    try:
        arguments.func(**dict((k, v) for k, v in vars(arguments).items()
            if k not in ('func', 'subcommand')))
    except ValueError as error:
        parser.error(error)


if __name__ == '__main__':
    main()
