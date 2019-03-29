#!/usr/bin/env python3
"""
@author: Timothy Baker

samplesheet_parser.py

"""

import pandas as pd


class SampleSheetParser:
    """

    A SampleSheetParser object that parses the input sample sheet in specific
    format and holds all instantiated paths.

    TODO:
        - implement logging
        - implement more secure error handling especially if sample sheet changes
        - devise better parsing methods
        - configure how this object talks with the docker container

    params:
        sample_sheet : ab/path/to/sample_sheet.csv

    attributes:
        path_info : dict of paths fastq_r1, fastq_r2, ref.fa, transcript.gtf
        adapters: dict of 5' and 3' adapter sequences for trimming
        barcode_seq : pandas dataframe of experiment design information and barcodes

    methods:
        parse_sample_sheet() : parses sample sheet and instantiates the attributes
        create_adapter_whitelist() : creates the barcode white list for zUMI from barcode_seq dict
        return_path_info() : returns path_info attr dict
        return_adapters() : return adapters dict
        return_barcode_seq() : returns pandas dataframe


    """

    def __init__(self, sample_sheet):
        self.sample_sheet = sample_sheet
        self.path_info = {'trimmed_r1': None, 'trimmed_r2': None}
        self.adapter = {}
        self.cell_data = None


    def parse_sample_sheet(self):
        """ parsing sample sheet and instantiating attributes """

        with open(self.sample_sheet, 'r') as csv_handle:

            line = csv_handle.readline()

            while line:

                # grabbing HEADER position
                if line.startswith('[HEADER]'):
                    header_offset = csv_handle.tell()

                # grabbing setting position
                if line.startswith('[SETTINGS]'):
                    settings_offset = csv_handle.tell()

                # grabbing data position
                if line.startswith('[DATA]'):
                    data_offset = csv_handle.tell()

                line = csv_handle.readline()

            # hacked way of parsing this csv, need to handle better
            # create header byte load
            byte_load = (settings_offset - header_offset) - 28

            # change position to header offset
            csv_handle.seek(header_offset)

            for line in csv_handle.readlines(byte_load):
                line_lst = line.split(',')

                # make all strings lower for slight error handling
                if line_lst[0].lower() == 'fastq_read1':
                    self.path_info['fastq_read1'] = line_lst[1]

                if line_lst[0].lower() == 'fastq_read2':
                    self.path_info['fastq_read2'] = line_lst[1]

                if line_lst[0].lower() == 'ref_genome':
                    self.path_info['ref_genome'] = line_lst[1]

                if line_lst[0].lower() == 'annotation':
                    self.path_info['annotation'] = line_lst[1]

                if line_lst[0].lower() == 'basename':
                    self.path_info['basename'] = line_lst[1]

            # setting trimmed fastq file paths
            # may need to fix paths
            # trimmed fastq goes to current directory, will need to find it
            # may not be able to use relative paths in docker
            self.path_info['trimmed_r1'] = './{}.trimmed.R1.fastq.gz'.format(\
                                                    self.path_info['basename'])

            self.path_info['trimmed_r2'] = './{}.trimmed.R2.fastq.gz'.format(\
                                                    self.path_info['basename'])

            # changing file position to settings/adapter offset
            csv_handle.seek(settings_offset)

            # setting adapter attributes
            adapter_list = csv_handle.readline().split(',')
            self.adapter['adapter_3'] = adapter_list[1]
            self.adapter['adapter_5'] = adapter_list[2]

            # changing file position to data offset to load into dataframe
            csv_handle.seek(data_offset)

            # pandas dataframe for easy
            self.cell_data = pd.read_csv(csv_handle)


    def create_adapter_whitelist(self):
        """ create barcode_whitelist text file """

        # need to direct toward specific path
        # will only exist in Docker
        self.cell_data.to_csv('barcode_white.txt', \
                                sep='\n', \
                                columns=['barcode_sequence'], \
                                header=False, \
                                index=False)

    def return_path_info(self):
        """ return path_info dict """
        return self.path_info

    def return_adapters(self):
        """ return adapter dict sequences """
        return self.adapter

    def return_barcode_seq(self):
        """ return cell barcode list """
        return self.cell_data



def main():
    """ run main for testing """
    pass

if __name__ == '__main__':
    main()
