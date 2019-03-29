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
from samplesheet_parser import SampleSheetParser


LOG_DIR = 'logs'
INITIAL_FASTQC_OUTPUT = 'logs/initial_fastqc'
TRIMMED_FASTQC_OUTPUT = 'logs/trimmed_fastqc'


if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)
    print("Made logs directory.")

if not os.path.exists(INITIAL_FASTQC_OUTPUT):
    os.makedirs(INITIAL_FASTQC_OUTPUT)
    print("Made initial fastqc output directory.")

if not os.path.exists(TRIMMED_FASTQC_OUTPUT):
    os.makedirs(TRIMMED_FASTQC_OUTPUT)
    print("Made trimmed fastqc output directory.")


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
    parser.add_argument('sample_sheet', type=str, help='Enter absolute/path/to/sample_sheet.csv')

    return parser.parse_args()



def build_fastqc_args(*args):
    """ takes fastqc args and creates a kwarg for fastqc command line.
        this is a nonexhaustive list of args. Would need to change to include
        new args.

        args must be in order of fastqc_keyword
    """

    # error handling
    assert isinstance(args[0], int), "Thread argument requires an int object."
    assert isinstance(args[1], str), "Output dir argument requires a str object."
    assert isinstance(args[2], str), "R1 input file argument requires a str object."
    assert isinstance(args[3], str), "R2 input file argument requires a str object."

    fastqc_keywords = ('threads', \
                        'output_dir', \
                        'fastq_read1', \
                        'fastq_read2'\
    )

    return dict(zip(fastqc_keywords, args))



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
    print("Running FASTQC.")
    fastqc_process = subprocess.Popen(fastqc_formatted_args, stdout=subprocess.PIPE, \
                                                        stderr=subprocess.PIPE)
    data_out, data_err = fastqc_process.communicate()

    # may need to fix output log location
    with open('./logs/fastqc-stdout.txt', 'wb') as output:
        output.write(data_out)
        output.write(data_err)



def build_cutadapt_args(*args):
    """ builds arg dict for running cutadapt; nonexhaustive list of args
        passable. will need to change this to include more if desired.

        args must be in in order of the cutadapt keywords
    """

    # error handling
    assert isinstance(args[0], int), "Thread argument requires an int object."
    assert isinstance(args[1], str), "3' Adapter argument requires a str object."
    assert isinstance(args[2], str), "5' Adapter argument requires a str object."
    assert isinstance(args[3], int), "Min length R1 argument requires an int object."
    assert isinstance(args[4], int), "Min length R2 argument requires an int object."
    assert isinstance(args[5], str), "R1 output file argument requires a str object."
    assert isinstance(args[6], str), "R2 output file argument requires a str object."
    assert isinstance(args[7], str), "R1 input file argument requires a str object."
    assert isinstance(args[8], str), "R2 input file argument requires a str object."

    # these args are not inclusive of all passable args
    # would need to fix this if other args are wanted
    cutadapt_keywords = ('threads', \
                        '3p_adapter', \
                        '5p_adapter', \
                        'min_length_r1', \
                        'min_length_r2', \
                        'output_read1', \
                        'output_read2', \
                        'fastq_read1', \
                        'fastq_read2' \
    )

    return dict(zip(cutadapt_keywords, args))



def run_cutadapt(**kwargs):
    """ runs cutadapt with provided arguments from build_cutadapt_args()
        params:
            kwargs : dict specified in build_cutadapt_args
        returns:
            None
    """

    cutadapt_command = """cutadapt
                        -j {threads}
                        -a {3p_adapter}
                        -A {5p_adapter}
                        -m {min_length_r1}:{min_length_r2}
                        -o {output_read1}
                        -p {output_read2}
                        {fastq_read1}
                        {fastq_read2}""".format(**kwargs)


    cutadapt_formatted_args = shlex.split(cutadapt_command)

    print("Running cutadapt.")

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

    LOGGER.info("Input args: %s", args)

    sample_sheet = args.sample_sheet
    threads = args.threads

    LOGGER.info("Created SampleSheetParser Object")
    sample_obj = SampleSheetParser(sample_sheet)

    sample_obj.parse_sample_sheet()
    LOGGER.info("Parsed sample sheet.")

    # creating barcode white list text file for zUMI
    sample_obj.create_adapter_whitelist()
    LOGGER.info("Created Barcode whitelist for zUMI.")

    # dict contains all relevant file paths
    file_path_info = sample_obj.return_path_info()

    # dict contains adapter trimming sequences
    adapter_info = sample_obj.return_path_info()


    initial_fastqc_kwarg = build_fastqc_args(threads, \
                                    INITIAL_FASTQC_OUTPUT, \
                                    file_path_info['fastq_read1'], \
                                    file_path_info['fastq_read2'])
    LOGGER.info("Building initial fastqc dictionary: %s", initial_fastqc_kwarg)
    LOGGER.info("Running FASTQC on initial raw data.")
    run_fastqc(**initial_fastqc_kwarg)

    cutadapt_kwargs = build_cutadapt_args(threads, \
                                        adapter_info['adapter_3'], \
                                        adapter_info['adapter_5'], \
                                        12, \
                                        20, \
                                        file_path_info['trimmed_r1'], \
                                        file_path_info['trimmed_r2'], \
                                        file_path_info['fastq_read1'], \
                                        file_path_info['fastq_read2'])
    LOGGER.info("Building cutadapt dictionary: %s", cutadapt_kwargs)
    LOGGER.info("Running cutadapt.")
    run_cutadapt(**cutadapt_kwargs)

    trimmed_fastqc_kwarg = build_fastqc_args(threads, \
                                    TRIMMED_FASTQC_OUTPUT, \
                                    file_path_info['trimmed_r1'], \
                                    file_path_info['trimmed_r2'])
    LOGGER.info("Building trimmed fastqc dictionary: %s", trimmed_fastqc_kwarg)
    LOGGER.info("Running FASTQC on trimmed data.")
    run_fastqc(**trimmed_fastqc_kwarg)



if __name__ == '__main__':
    main()
