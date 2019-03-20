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
import pathlib
import argparse
import logging
import numpy as np
import pandas as pd
from Bio import SeqIO
import time



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
    parser.add_argument('sample_sheet', type=str, help='Enter absolute/path/to/sample_sheet.csv')
    parser.add_argument('fastq_r1', type=str, help='Enter absolute/path/to/fastq_read1')
    parser.add_argument('fastq_r2', type=str, help='Enter absolute/path/to/fastq_read2')

    return parser.parse_args()




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



def handle_sample_sheet(sample_sheet):
    """ takes sample sheet and returns a dict based on grouping. """

    try:
        # create new directory that stores demultiplexed fastqs
        os.mkdir('fastq_demultiplex')
    except FileExistsError:
        pass

    # load sample sheet into pandas dataframe
    sample_sheet_df = pd.read_csv(sample_sheet)

    # instantiate new columns with default values
    sample_sheet_df.assign(**{'barcode_count': 0, \
                            'barcode_hash': None, \
                            'umi_sequence': None, \
                            'fastq_r1_path': None, \
                            'fastq_r2_path': None})

    try:
        for index, row in sample_sheet_df.iterrows():

            # creating hash values for barcode comparison
            sample_sheet_df.at[index, 'barcode_hash'] = hash(row['barcode_sequence'])

            # initializing demultiplexed fastq files based on specified experiment group
            fq_r1_name = './fastq_demultiplex/fastq_R1.group_{}.demultiplexed.fastq'.format(row['experiment_group'])
            fq_r2_name = './fastq_demultiplex/fastq_R2.group_{}.demultiplexed.fastq'.format(row['experiment_group'])
            fq_r1_und_name = './fastq_demultiplex/fastq_R1.group_{}.undetermined.demultiplexed.fastq'.format(row['experiment_group'])
            fq_r2_und_name = './fastq_demultiplex/fastq_R2.group_{}.undetermined.demultiplexed.fastq'.format(row['experiment_group'])

            # creating empty files to append to later
            open(fq_r1_name, 'a+').close()
            open(fq_r2_name, 'a+').close()
            open(fq_r1_und_name, 'a+').close()
            open(fq_r2_und_name, 'a+').close()

            # using Path object to grab absolute path to store in df
            fq_r1 = str(pathlib.Path(fq_r1_name).resolve())
            fq_r2 = str(pathlib.Path(fq_r2_name).resolve())
            fq_r1_und = str(pathlib.Path(fq_r1_und_name).resolve())
            fq_r2_und = str(pathlib.Path(fq_r2_und_name).resolve())

            sample_sheet_df.at[index, 'fastq_r1_path'] = fq_r1
            sample_sheet_df.at[index, 'fastq_r2_path'] = fq_r2
            sample_sheet_df.at[index, 'fastq_r1_und_path'] = fq_r1_und
            sample_sheet_df.at[index, 'fastq_r2_und_path'] = fq_r2_und

    except OSError as error:
        # need to add more error handling
        print(error)


    return sample_sheet_df




def demultiplex_fastq_files(fastq_r1, fastq_r2, sample_sheet_df):


    # loads path of index db
    fastq_r1_idx_db = './fastq_demultiplex/fastq_r1_idx_db.idx'
    fastq_r2_idx_db = './fastq_demultiplex/fastq_r2_idx_db.idx'

    # creates a sqllite3 index db
    fastq_r1_record_db = SeqIO.index_db(fastq_r1_idx_db, fastq_r1, 'fastq')
    fastq_r2_record_db = SeqIO.index_db(fastq_r2_idx_db, fastq_r2, 'fastq')


    for record in fastq_r1_record_db:

        # use get_raw() for access to byte-like object
        fastq_r1_seq = fastq_r1_record_db.get_raw(record).decode().split('\n')[1]

        # parsing out UMI first 6 base pairs of read 1 sequence
        sample_sheet_df['umi_sequence'] = fastq_r1_seq[0:6]

        # parsing cel-barcode and creating a hash for fast comparison
        # cel-barcode is 6-mer
        fastq_r1_seq_hash = hash(fastq_r1_seq[6:12])

        match_row = sample_sheet_df.loc[\
                    sample_sheet_df.barcode_hash \
                    == fastq_r1_seq_hash].to_numpy().flatten()

        mismatch_row = sample_sheet_df.loc[\
                        sample_sheet_df.barcode_hash \
                        != fastq_r1_seq_hash].to_numpy().flatten()

        fastq_r1_raw_byte = fastq_r1_record_db.get_raw(record)
        fastq_r2_raw_byte = fastq_r2_record_db.get_raw(record)

        if match_row.any():
            r1_file = open(str(match_row[6]), 'ab')
            r2_file = open(str(match_row[7]), 'ab')
            r1_file.write(fastq_r1_raw_byte)
            r2_file.write(fastq_r2_raw_byte)
        else:
            r1_und_file = open(str(mismatch_row[8]), 'ab')
            r2_und_file = open(str(mismatch_row[9]), 'ab')
            r1_und_file.write(fastq_r1_raw_byte)
            r2_und_file.write(fastq_r2_raw_byte)



def main():
    """ creates the tmp directories and splits the fastq into smaller
        files for increased iteration
    """

    args = arg_parser()

    sample_sheet = args.sample_sheet
    fastq_r1 = args.fastq_r1
    fastq_r2 = args.fastq_r2

    sample_sheet_df = handle_sample_sheet(sample_sheet)

    demultiplex_fastq_files(fastq_r1, fastq_r2, sample_sheet_df)

if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print(end - start)
