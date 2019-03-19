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
import logging
import HTSeq

import pandas as pd


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
