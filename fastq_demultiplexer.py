#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

fastq_demultiplexer.py

# grabs the directory the python file is in
# dir_path = os.path.dirname(os.path.realpath(__file__))
# need to grab and open fastqc unzipped output
# create script that does this
# need to handle gunzipped fastqs

# need to deal with zipped fastqs! test with unzipped!

"""

import os
import sys
import argparse
import logging
import collections
import HTSeq
import pandas as pd
import pprint as pp
from Bio import SeqIO


LOGGER = logging.getLogger(__name__)

LOGGER.setLevel(logging.INFO)

FORMATTER = logging.Formatter('%(levelname)s:%(name)s:%(asctime)s:%(message)s')

FILE_HANDLER = logging.FileHandler("fastq_demultiplex.log")

FILE_HANDLER.setFormatter(FORMATTER)

LOGGER.addHandler(FILE_HANDLER)


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Demultiplexes single cell RNA-seq FASTQs by group\n'
    )
    # add statistic file for logs
    # parser.add_argument('-o', '--out-dir', metavar='DIRNAME', type=str, help='path to text file that has ftp links')
    # parser.add_argument('sample_sheet', type=str, help='Enter absolute/path/to/sample_sheet.csv')
    parser.add_argument('fastq_r1', type=str, help='Enter absolute/path/to/fastq_read1')
    parser.add_argument('fastq_r2', type=str, help='Enter absolute/path/to/fastq_read2')

    return parser.parse_args()




 #   _____ _____  _      _____ _______   ______       _____ _______ ____
 #  / ____|  __ \| |    |_   _|__   __| |  ____/\    / ____|__   __/ __ \
 # | (___ | |__) | |      | |    | |    | |__ /  \  | (___    | | | |  | |
 #  \___ \|  ___/| |      | |    | |    |  __/ /\ \  \___ \   | | | |  | |
 #  ____) | |    | |____ _| |_   | |    | | / ____ \ ____) |  | | | |__| |
 # |_____/|_|    |______|_____|  |_|    |_|/_/    \_|_____/   |_|  \___\_\
 #

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def create_fastqtmp_dirs():
    """ grabs current working directory and creates a top level tmp directory
        and two subdirectories called fastq_r1, and fastq_r2.
    """

    cur_dir = os.getcwd()

    # ensuring absolute path for traceability
    fastq_r1_tmp_path = cur_dir + '/fastqtmp/fastq_r1'
    fastq_r2_tmp_path = cur_dir + '/fastqtmp/fastq_r2'

    # change the prints to logging
    try:
        os.makedirs(fastq_r1_tmp_path)
        os.makedirs(fastq_r2_tmp_path)
    except OSError:
        LOGGER.info("Failed to create %s and %s directories.", fastq_r1_tmp_path, \
                                                                fastq_r2_tmp_path)
    else:
        LOGGER.info("Successfully created the directory.")




def detect_fastq_file(fastq):
    """ takes fastq path and detects whether it is R1 or R2 """

    LOGGER.info("Path to fastq %s", str(fastq))

    # need to ensure that the path has no other R1 or R2 present
    if 'R1' in fastq.upper():
        return 1
    elif 'R2' in fastq.upper():
        return 2
    return False
        # log to file if this occurs



def split_fastq_files(fastq_file, total_seq_num):
    """ takes fastq and splits them up by total_record_num / 25
        cores.
        fastq (str) : absolute path to fastq file
        total_seq_num (int) : number of sequences from fastqc report
    """

    # loads the fastq into a SeqIO object
    fastq_record = SeqIO.parse(open(fastq_file), "fastq")

    # add assert statements
    # may need to consider other machines that do not have enough cores
    num_of_records_to_split = int(total_seq_num) / 25

    # this has to be fastq path, need to ensure path has no other
    # R1 or R2 present
    fastq_num = detect_fastq_file(fastq_file)
    LOGGER.info("FASTQ R%s detected", str(fastq_num))

    for idx, batch in enumerate(batch_iterator(fastq_record, num_of_records_to_split)):

        cur_dir = os.getcwd()

        fastq_r1_tmp_path = cur_dir + '/fastqtmp/fastq_r1/'

        fastq_filename = "multiplex.fastq_R{}_{}.fastq".format(fastq_num, idx + 1)

        abs_filename = fastq_r1_tmp_path + fastq_filename

        with open(abs_filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fastq")

        # change to logging function
        print("Created {}".format(abs_filename))
        LOGGER.info("Wrote %s records to %s", str(count), abs_filename)



def handle_sample_sheet(sample_sheet):
    """ takes sample sheet and returns a dict based on grouping. """

    # Need to add a lot of checks to ensure data is correct
    # may need to convert to dict rather than dataframe

    exp_group_dict = {}

    sample_sheet_df = pd.read_csv(sample_sheet)

    experiment_group_num = set([int(x) for x in sample_sheet_df['experiment_group']])

    sample_sheet_group = sample_sheet_df.groupby(['experiment_group'])

    # how would I access if i changed from dataframe to dict, probably better on memory
    for condition in experiment_group_num:
        exp_group_dict[condition] = sample_sheet_group.get_group(condition)

    return exp_group_dict


# class FastqDemultiplexer:
#
#     def __init__(self, barcode_index_file, sample_sheet, fastq_r1, fastq_r2):
#         self.sample_sheet_file = sample_sheet
#         self.fastq_r1 = fastq_r1
#         self.fastq_r2 = fastq_r2
#
#
#     def build_barcode_dict(self):
#
#         barcode_dict = {}
#
#         with open(self.barcode_index_file, 'rb') as bar_idx:
#



def main():
    """ creates the tmp directories and splits the fastq into smaller
        files for increased iteration
    """

    args = arg_parser()

    fastq_r1 = args.fastq_r1
    fastq_r2 = args.fastq_r2

    create_fastqtmp_dirs()

    total_seq_num = 377411561

    split_fastq_files(fastq_r1, total_seq_num)

    split_fastq_files(fastq_r2, total_seq_num)




if __name__ == '__main__':
    main()
