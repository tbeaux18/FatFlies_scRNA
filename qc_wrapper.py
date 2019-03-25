#!/usr/bin/env python3
"""
@author: Timothy Baker

wrapper_functions.py

"""

import os
import argparse
import shlex
import subprocess


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Python3 wrapper for quality control on single cell RNA-seq'
    )
    parser.add_argument('-a', '--adapter_file', type=str, help='absolute path to the adapter file')
    parser.add_argument('-c', '--contam_file', type=str, help='absolute path to the contamination file')
    parser.add_argument('-b', '--basename', type=str, help='basename to trimmed fastq files')
    parser.add_argument('fastq_r1', type=str, help='Enter absolute/path/to/fastq_read1')
    parser.add_argument('fastq_r2', type=str, help='Enter absolute/path/to/fastq_read2')

    return parser.parse_args()



def build_fastqc_arg(*args):
    """ Takes x number of args and formats to fastqc_arg for downstream input.
        No more than 3 args, only necessary to properly format trim_galore.
        const_params:
            --extract : auto unzip output files
            -t/--threads : 4
            -k/--kmer : 6
        params:
            adapter (str) : ab/path/to/adapter_file
            contam (str) : ab/path/to/contam_file
            fastqc_out (str) : ab/path/to/fastqc_output
                                requires directory to be made
        returns:
            fastqc_arg (str) : full formatted str object for fastqc
    """

    # raises assertion error if not true
    assert len(args) == 3, "3 args were not passed."


    # leave kmer constant as 6 for cel-seq2
    # umis and barcodes are 6 bp long
    fastqc_arg = "--extract --threads 4 --kmers 6 --format fastq" \
                "--adapters {} --contaminants {} --outdir {}".format(*args)

    return fastqc_arg



def build_trim_dict(*args):
    """ Builds the trim_galore kwarg dict to pass into the trim command function.
        Only a max of 5 args should be passed, and order matters. Only pass absolute
        paths, not optimized for relative paths.
        Keeps fastqs zipped
        params:
            fastqc_arg (str) : fastqc_arg command to run fastqc
            basename (str) : basename of trimmed fastq out put
            trim_output (str) : /path/to/trim_output
            fastq_read1 (str) : /path/to/untrimmed_fastq_read1.fastq.gz
            fastq_read2 (str) : /path/to/untrimmed_fastq_read2.fastq.gz
        returns:
            trimgalore_dict (dict)
    """

    # raise error if not true
    assert len(args) == 5, "5 args were not passed."

    # constant keyword args
    trim_galore_keys = ('fastqc_arg', 'basename', 'trim_output', 'fastq_read1', 'fastq_read2')

    # zips the keyword and arg paths into a dict for downstream analysis
    # order is imperative
    trimgalore_dict = dict(zip(trim_galore_keys, args))

    return trimgalore_dict



def run_trim_galore(**kwargs):
    """ Takes a keyword dict and passes into the trimgalore command and calls
        a process to begin while piping standard output and standard error.
        constant_params:
            -f/--format : fastq
            --small_rna : will trim the following sequences from
                            R1: TGGAATTCTCGG & R2: GATCGTCGGACT
            -j/--cores : 4
            -s/--stringency : 2
        params:
            trimgalore_dict (dict) : dict with all keywords and values for
                    fastqc_args, basename, trim_output, fastq_read1, fastq_read2
        returns:
            None
    """

    trimgalore_command = """trim_galore --fastqc_args "{fastqc_arg}"
                        --basename {basename}
                        --small_rna
                        -j 4
                        --stringency 2
                        -o {trim_output}
                        --paired {fastq_read1} {fastq_read2}""".format(**kwargs)

    trimgalore_formatted_args = shlex.split(trimgalore_command)

    trimgalore_process = subprocess.Popen(trimgalore_formatted_args, \
                                            stdout=subprocess.PIPE, \
                                            stderr=subprocess.PIPE)
    data_out, data_err = trimgalore_process.communicate()

    with open('trimgalore-stdout.txt', 'ab') as output:
        output.write(data_out)
        output.write(data_err)




def main():
    """ initialize command line args, build trim_galore command, and execute """

    # initializing arg parser
    args = arg_parser()

    # setting command line arguments to be passed
    fastq_read1 = args.fastq_r1
    fastq_read2 = args.fastq_r2
    basename = args.basename
    adapter_file = args.adapter_file
    contam_file = args.contam_file


    # create the fastqc_output path
    # need to track directory structure
    fastqc_out_dir = 'trimgal_fastqc'
    trimgal_out_dir = 'trimgal_output'
    if not os.path.exists(fastqc_out_dir):
        os.mkdir(fastqc_out_dir)
        os.mkdir(trimgal_out_dir)


    fastqc_output = './trimgal_fastqc'
    trim_output = './trimgal_output'

    # order is imperative
    # 1) adapter_file 2) contam_file 3) fastqc_output path
    fastqc_out = build_fastqc_arg(adapter_file, contam_file, fastqc_output)


    # 1) fastqc_out path 2) basename 3) trim output dir 4) fastq R1 5) fastq R2
    trim_kwarg = build_trim_dict(fastqc_out, basename, trim_output, fastq_read1, fastq_read2)

    # executes the script
    # need to include error handling and log outputs
    run_trim_galore(**trim_kwarg)

if __name__ == '__main__':
    main()
