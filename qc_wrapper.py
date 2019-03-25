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
    parser.add_argument('fastq_r1', type=str, help='Enter absolute/path/to/fastq_read1')
    parser.add_argument('fastq_r2', type=str, help='Enter absolute/path/to/fastq_read2')

    return parser.parse_args()



def run_cutadapt(fastq_read1, fastq_read2, celseq=True):
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


    cutadapt_output = 'cutadapt_output'

    if not os.path.exists(cutadapt_output):
        os.mkdir(cutadapt_output)

    trimmed_read1_fastq = './cutadapt_output/cutadapt.trimmed.R1.fastq.gz'
    trimmed_read2_fastq = './cutadapt_output/cutadapt.trimmed.R1.fastq.gz'

    if celseq:
        r1_adapter_sequence = 'TGGAATTCTCGG'
        r2_adapter_sequence = 'GATCGTCGGACT'
    else:
        pass # fix later


    cutadapt_keywords = ('read1_adapter', 'read2_adapter', \
                        'read1_output_file', 'read2_output_file', \
                        'fastq_read1', 'fastq_read2')

    cutadapt_args = (r1_adapter_sequence, r2_adapter_sequence, \
                    trimmed_read1_fastq, trimmed_read2_fastq, \
                    fastq_read1, fastq_read2)

    cutadapt_kwargs = dict(zip(cutadapt_keywords, cutadapt_args))

    cutadapt_command = """cutadapt
                        -j 4
                        -a {read1_adapter}
                        -A {read2_adapter}
                        -m 12:20
                        -o {read1_output_file}
                        -p {read2_output_file}
                        {fastq_read1}
                        {fastq_read2}""".format(**cutadapt_kwargs)


    cutadapt_formatted_args = shlex.split(cutadapt_command)
    print("Running cutadapt.")

    cutadapt_process = subprocess.Popen(cutadapt_formatted_args, \
                                            stdout=subprocess.PIPE, \
                                            stderr=subprocess.PIPE)
    data_out, data_err = cutadapt_process.communicate()

    with open('cutadapt-stdout.txt', 'ab') as output:
        output.write(data_out)
        output.write(data_err)


    return trimmed_read1_fastq, trimmed_read2_fastq



def run_fastqc(*args):
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
    assert len(args) == 4, "4 args were not passed."

    fastqc_out_dir = 'fastqc_run'
    if not os.path.exists(fastqc_out_dir):
        os.mkdir(fastqc_out_dir)

    # leave kmer constant as 6 for cel-seq2
    # umis and barcodes are 6 bp long
    # space at the end of the first string is crucial to fastqc operating
    fastqc_command = "fastqc --extract --threads 4 --kmers 6 --format fastq " \
                "--adapters {} --contaminants {} --outdir ./fastqc_run {} {}".format(*args)

    fastqc_formatted_args = shlex.split(fastqc_command)
    print("Running FASTQC.")
    fastqc_process = subprocess.Popen(fastqc_formatted_args, stdout=subprocess.PIPE, \
                                                        stderr=subprocess.PIPE)
    data_out, data_err = fastqc_process.communicate()

    with open('fastqc-stdout.txt', 'wb') as output:
        output.write(data_out)
        output.write(data_err)


def main():
    """ initialize command line args, build trim_galore command, and execute """

    # initializing arg parser
    args = arg_parser()

    # setting command line arguments to be passed
    fastq_read1 = args.fastq_r1
    fastq_read2 = args.fastq_r2
    adapter_file = args.adapter_file
    contam_file = args.contam_file

    # executes the script
    # need to include error handling and log outputs
    r1_trim_fastq, r2_trim_fastq = run_cutadapt(fastq_read1, fastq_read2)

    # order is imperative
    # 1) adapter_file 2) contam_file 3) read1 trimmed path 4) read2 trimmed path
    run_fastqc(adapter_file, contam_file, r1_trim_fastq, r2_trim_fastq)

if __name__ == '__main__':
    main()
