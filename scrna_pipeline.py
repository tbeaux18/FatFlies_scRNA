#!/usr/bin/env python3
"""
@authors: Timothy Baker


scrna_pipeline.py

script runs the main pipeline

"""

import os
import pathlib
import csv
import shlex
import argparse
import subprocess
import samplesheet_parser


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Runs the single-cell RNA-seq pipeline.'
    )
    parser.add_argument('sample_sheet', type=str, help='Enter absolute/path/to/sample_sheet.csv')

    return parser.parse_args()


def run_qc(fastq_read1, fastq_read2):
    """ takes fastq files and starts the qc wrapper """

    run_qc_wrapper = "python3 qc_wrapper.py {} {}".format(fastq_read1, fastq_read2)

    run_qc_formatted_args = shlex.split(run_qc_wrapper)

    qc_wrapper_process = subprocess.Popen(run_qc_formatted_args, stdout=subprocess.PIPE, \
                                                            stderr=subprocess.PIPE)
    data_out, data_err = qc_wrapper_process.communicate()

    print(data_out, data_err)



def run_zumi(zumi_config):
    """ takes zumi config file and initates zumi pipeline """

    # need to track zumi directory
    zumi_command = "zumi-master.sh -y {}".format(zumi_config)

    zumi_formatted_args = shlex.split(zumi_command)

    zumi_process = subprocess.Popen(zumi_formatted_args, stdout=subprocess.PIPE, \
                                                            stderr=subprocess.PIPE)
    data_out, data_err = zumi_process.communicate()

    with open('zumi-stdout.txt', 'wb+') as zumi_output:
        zumi_output.write(data_out)
        zumi_output.write(data_err)

def main():
    """ runs main """

    args = arg_parser()

    sample_sheet = args.sample_sheet

    # instantiated object
    sample_sheet_object = samplesheet_parser.SampleSheetParser(sample_sheet)

    sample_sheet_object.parse_sample_sheet()
    fastq_read1 = sample_sheet_object.return_read1()
    fastq_read2 = sample_sheet_object.return_read2()

    # run the qc wrapper.py
    run_qc(fastq_read1, fastq_read2)



if __name__ == '__main__':
    main()
