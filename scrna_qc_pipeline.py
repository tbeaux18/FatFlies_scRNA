#!/usr/bin/env python3
"""
@authors: Timothy Baker
@version: 1.0.0

scrna_pipeline_qc.py

Script takes in sample sheet and thread count, creates a SampleSheetParser object
which stores all information from the sample sheet. Using the information,
quality control is performed on the raw data, and trimmed data. There has been no
implementation of error handling, or parsing the output log files before
continuing onto the next step.

TODO:
    Parsing log files
    Error handling
    Path name control

Dependencies:
    python3
    samplesheet_parser.py must be in same root directory
    fastqc
    cutadapt -- installed for python3 due to multithreading

"""

import os
import shlex
import argparse
import subprocess
import logging



INITIAL_FASTQC_OUTPUT = 'logs/initial_fastqc'
TRIMMED_FASTQC_OUTPUT = 'logs/trimmed_fastqc'

if not os.path.exists(INITIAL_FASTQC_OUTPUT):
    os.makedirs(INITIAL_FASTQC_OUTPUT)


if not os.path.exists(TRIMMED_FASTQC_OUTPUT):
    os.makedirs(TRIMMED_FASTQC_OUTPUT)



LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
FORMATTER = logging.Formatter('%(levelname)s:%(name)s:%(asctime)s:%(message)s')
FILE_HANDLER = logging.FileHandler("logs/qc_log.log")
FILE_HANDLER.setFormatter(FORMATTER)
LOGGER.addHandler(FILE_HANDLER)



def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Runs the single-cell RNA-seq pipeline.'
    )
    parser.add_argument('-t', '--threads', type=int, help='Enter thread int object')

    parser.add_argument('-f', '--fastq_r1', type=str, help='Enter abs/path/to/fastq_r1')
    parser.add_argument('-p', '--fastq_r2', type=str, help='Enter abs/path/to/fastq_r2')

    parser.add_argument('-i', '--trimmed_r1', type=str, help='Enter abs/path/to/trimmed_r1 output')
    parser.add_argument('-d', '--trimmed_r2', type=str, help='Enter abs/path/to/trimmed_r2 output')

    parser.add_argument('-a', '--adapter_3', type=str, help='Enter 3p adapter sequence')
    parser.add_argument('-A', '--adapter_5', type=str, help='Enter 5p adapter sequence')

    parser.add_argument('-m', '--read1_min', type=int, help='Enter int for min length for read 1')
    parser.add_argument('-M', '--read2_min', type=int, help='Enter int for min length for read 2')

    return parser.parse_args()




def parse_fastqc_results(fastqc_path):
    """ parses the fastqc summary file and stores results

        There are two sub directories to parse from, the initial
        and the trimmed.

        Initial:
            ./logs/initial_fastqc/<basename>_fastqc/summary.txt
        Trimmed:
            ./logs/trimmed_fastqc/<basename>_fastqc/summary.txt
        params:
            fastqc_path : path to summary.txt
        returns:
            summary_results : dict of results
    """

    summary_results = {}

    # needs to be abs/path/to/summary.txt
    with open(fastqc_path, 'r') as fastqc_file:
        for line in fastqc_file:
            line = line.strip().split('\t')
            summary_results[line[1]] = line[0]

    return summary_results



def run_fastqc(**kwargs):
    """ runs fastqc with the provided args from build_fastqc_args()
        params:
            kwargs : dict specified in build_fastqc_args
            --extract : qc log files will be unzipped
            --format : fastq default
        returns:
            None
    """

    fastqc_command = """fastqc
                        --extract
                        --threads {threads}
                        --format fastq
                        --outdir {output_dir}
                        {fastq_read1}
                        {fastq_read2}""".format(**kwargs)

    fastqc_formatted_args = shlex.split(fastqc_command)

    LOGGER.info("FASTQC args: %s", fastqc_formatted_args)

    fastqc_process = subprocess.Popen(fastqc_formatted_args, stdout=subprocess.PIPE, \
                                                        stderr=subprocess.PIPE)
    data_out, data_err = fastqc_process.communicate()

    # may need to fix output log location
    with open('./logs/fastqc-stdout.txt', 'wb') as output:
        output.write(data_out)
        output.write(data_err)



def run_cutadapt(**kwargs):
    """ runs cutadapt with provided arguments from build_cutadapt_args()
        params:
            kwargs : dict specified in build_cutadapt_args
        returns:
            None
    """

    cutadapt_command = """cutadapt
                        -j {threads}
                        -a {adapter_3}
                        -A {adapter_5}
                        -m {min_length_r1}:{min_length_r2}
                        -o {output_read1}
                        -p {output_read2}
                        {fastq_read1}
                        {fastq_read2}""".format(**kwargs)


    cutadapt_formatted_args = shlex.split(cutadapt_command)

    LOGGER.info("cutadapt args: %s", cutadapt_formatted_args)

    cutadapt_process = subprocess.Popen(cutadapt_formatted_args, \
                                            stdout=subprocess.PIPE, \
                                            stderr=subprocess.PIPE)
    data_out, data_err = cutadapt_process.communicate()

    # may need to fix log location
    with open('./logs/cutadapt-stdout.txt', 'ab') as output:
        output.write(data_out)
        output.write(data_err)



def main():
    """ takes command line args, creates SampleSheetParser object, and runs
        quality control
    """

    args = arg_parser()

    LOGGER.info("Args passed to script: %s", args)

    threads = args.threads
    fastq_r1 = args.fastq_r1
    fastq_r2 = args.fastq_r2
    trimmed_r1 = args.trimmed_r1
    trimmed_r2 = args.trimmed_r2
    adapter_3 = args.adapter_3
    adapter_5 = args.adapter_5
    read1_min = args.read1_min # 12
    read2_min = args.read2_min # 20


    LOGGER.info("Running FASTQC on initial raw data.")

    run_fastqc(
        threads=threads,
        output_dir=INITIAL_FASTQC_OUTPUT,
        fastq_read1=fastq_r1,
        fastq_read2=fastq_r2
    )




    LOGGER.info("Running cutadapt.")

    run_cutadapt(
        threads=threads,
        adapter_3=adapter_3,
        adapter_5=adapter_5,
        min_length_r1=read1_min,
        min_length_r2=read2_min,
        output_read1=trimmed_r1,
        output_read2=trimmed_r2,
        fastq_read1=fastq_r1,
        fastq_read2=fastq_r2
    )

    LOGGER.info("Running FASTQC on trimmed data.")

    run_fastqc(
        threads=threads,
        output_dir=TRIMMED_FASTQC_OUTPUT,
        fastq_read1=trimmed_r1,
        fastq_read2=trimmed_r2
    )



if __name__ == '__main__':
    main()
