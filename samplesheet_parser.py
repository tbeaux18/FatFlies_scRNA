#!/usr/bin/env python3
"""
@author: Timothy Baker

samplesheet_parser.py

04/01/2019:
        - Broke out original parse_sample_sheet() method to smaller
            parsing functions.
        - included new attributes for downstream configuration for zUMI
            and differential expression inputs
        - need to handle parsing a lot better; perhaps break into subclasses
"""

import pandas as pd


class SampleSheetParser:
    """

    A SampleSheetParser object that parses the input sample sheet in specific
    format and holds all instantiated paths.

    TODO:
        - fix diff_exp input, can only accept concatenated numbers
            and turns them into a list of strings
        - implement logging
        - implement more secure error handling especially if sample sheet changes
        - devise better parsing methods
        - configure how this object talks with the docker container

    params:
        sample_sheet : ab/path/to/sample_sheet.csv

    attributes:
        offset_pos : dict of offset positions for better troubleshooting of each [VALUE]
        header_info : dict of header info date, run_name, library_prep, basename
        path_info : dict of paths fastq_r1, fastq_r2, ref.fa, transcript.gtf
        zumi_input : dict of zumi params that will be added to config file, include
                        bc_filter_num_bases, bc_filter_phred, bc_ham_dist
                        umi_filter_num_bases, umi_filter_phred, umi_ham_dist
                        zum_start_stage
        diff_input : dict of condition groups to use during differential expression
        adapters: dict of 5' and 3' adapter sequences for trimming
        cell_data : pandas dataframe of experiment design information and barcodes

    methods:
        locate_offsets() : finds each [VALUE] offset position
        parse_sample_sheet_header() : stores header info
        parse_sample_sheet_zumi() : stores zumi input config info
        parse_sample_sheet_diffexp() : stores info for diff expression
        parse_sample_sheet_adapters() : stores adapter info for trimming
        parse_sample_sheet_data() : stores cell data
        run_parsing_methods() : runs all instance methods
        create_adapter_whitelist() : creates the barcode white list for zUMI from barcode_seq dict
        return_offsets() : returns offset position dict
        return_header_info() : returns header info
        return_zumi_input() : returns zumi config input info
        return_diff_input() : returns diff expression input info
        return_path_info() : returns path_info attr dict
        return_adapters() : return adapters dict for trimming
        return_cell_data() : returns pandas dataframe of cell info

    """

    def __init__(self, sample_sheet):
        self.sample_sheet = sample_sheet
        self.offset_pos = {}
        self.header_info = {}
        self.path_info = {}
        self.zumi_input = {}
        self.diff_input = {}
        self.adapter = {}
        self.cell_data = None


    def locate_offsets(self):
        """ parse sample sheet and grab offset positions of bracketed names """
        with open(self.sample_sheet, 'r') as csv_handle:

            line = csv_handle.readline()

            while line:

                # grabbing header position
                if line.startswith('[HEADER]'):
                    self.offset_pos['header_offset'] = csv_handle.tell()

                # grabbing zumi position
                if line.startswith('[ZUMI]'):
                    self.offset_pos['zumi_offset'] = csv_handle.tell()

                # grabbing diff_exp position
                if line.startswith('[DIFF_EXP]'):
                    self.offset_pos['diff_offset'] = csv_handle.tell()

                # grabbing setting position
                if line.startswith('[ADAPTERS]'):
                    self.offset_pos['adapters_offset'] = csv_handle.tell()

                # grabbing data position
                if line.startswith('[DATA]'):
                    self.offset_pos['data_offset'] = csv_handle.tell()

                line = csv_handle.readline()


    def parse_sample_sheet_header(self):
        """ parsing the header section of the sample sheet """

        with open(self.sample_sheet, 'r') as csv_handle:

            header_byte_load = (\
            self.offset_pos['zumi_offset'] - self.offset_pos['header_offset']\
            ) - 28

            csv_handle.seek(self.offset_pos['header_offset'])

            for line in csv_handle.readlines(header_byte_load):
                line_lst = line.split(',')

                # parsing overall project information to be used for logging, &
                # filenaming
                if line_lst[0].lower() == 'date':
                    self.header_info['date'] = line_lst[1]

                if line_lst[0].lower() == 'run_name':
                    self.header_info['run_name'] = line_lst[1]

                if line_lst[0].lower() == 'library_prep':
                    self.header_info['library_prep'] = line_lst[1]

                if line_lst[0].lower() == 'basename':
                    self.header_info['basename'] = line_lst[1]

                # parsing path information from header
                # grouping relevant information together
                if line_lst[0].lower() == 'fastq_read1':
                    self.path_info['fastq_read1'] = line_lst[1]

                if line_lst[0].lower() == 'fastq_read2':
                    self.path_info['fastq_read2'] = line_lst[1]

                if line_lst[0].lower() == 'ref_genome':
                    self.path_info['ref_genome'] = line_lst[1]

                if line_lst[0].lower() == 'annotation':
                    self.path_info['annotation'] = line_lst[1]


            # setting trimmed fastq file paths
            # may need to fix paths
            # trimmed fastq goes to current directory, will need to find it
            # may not be able to use relative paths in docker
            self.path_info['trimmed_r1'] = './{}.trimmed.R1.fastq.gz'.format(\
                                                    self.header_info['basename'])

            self.path_info['trimmed_r2'] = './{}.trimmed.R2.fastq.gz'.format(\
                                                    self.header_info['basename'])


    def parse_sample_sheet_zumi(self):
        """ parsing the zumi section of the sample sheet """

        with open(self.sample_sheet, 'r') as csv_handle:

            zumi_byte_load = (\
            self.offset_pos['diff_offset'] - self.offset_pos['zumi_offset']\
            ) - 28

            # changing file position to zumi offset
            csv_handle.seek(self.offset_pos['zumi_offset'])

            # setting up zumi dict for yaml config file generation
            for line in csv_handle.readlines(zumi_byte_load):
                line_lst = line.split(',')

                if line_lst[0].lower() == 'bc_filter_num_bases':
                    self.zumi_input['bc_filter_num_bases'] = line_lst[1]

                if line_lst[0].lower() == 'bc_filter_phred':
                    self.zumi_input['bc_filter_phred'] = line_lst[1]

                if line_lst[0].lower() == 'bc_ham_dist':
                    self.zumi_input['bc_ham_dist'] = line_lst[1]

                if line_lst[0].lower() == 'umi_filter_num_bases':
                    self.zumi_input['umi_filter_num_bases'] = line_lst[1]

                if line_lst[0].lower() == 'umi_filter_phred':
                    self.zumi_input['umi_filter_phred'] = line_lst[1]

                if line_lst[0].lower() == 'umi_ham_dist':
                    self.zumi_input['umi_ham_dist'] = line_lst[1]

                if line_lst[0].lower() == 'zum_start_stage':
                    self.zumi_input['zum_start_stage'] = line_lst[1]


    def parse_sample_sheet_diffexp(self):
        """ parsing the diff_expression section of the sample sheet """

        with open(self.sample_sheet, 'r') as csv_handle:

            diffexp_byte_load = (\
            self.offset_pos['adapters_offset'] - self.offset_pos['diff_offset']\
            ) - 28

            csv_handle.seek(self.offset_pos['diff_offset'])

            for line in csv_handle.readlines(diffexp_byte_load):
                line_lst = line.split(',')
                if line_lst[0].lower() == 'test_group':
                    # need to fix this
                    self.diff_input['test_group'] = list(line_lst[1])

                if line_lst[0].lower() == 'control_group':
                    self.diff_input['control_group'] = list(line_lst[1])


    def parse_sample_sheet_adapters(self):
        """ parsing the adapter section of sample sheet """

        with open(self.sample_sheet, 'r') as csv_handle:

            # changing to adapter offset position
            csv_handle.seek(self.offset_pos['adapters_offset'])

            # setting adapter attributes
            adapter_list = csv_handle.readline().split(',')

            # adapter sequences to trim
            self.adapter['adapter_3'] = adapter_list[1]
            self.adapter['adapter_5'] = adapter_list[2]


    def parse_sample_sheet_data(self):
        """ parsing cell data section from sample sheet """

        with open(self.sample_sheet, 'r') as csv_handle:

            # changing to data offset position
            csv_handle.seek(self.offset_pos['data_offset'])

            # pandas dataframe for easy storage, and retrieval
            self.cell_data = pd.read_csv(csv_handle)


    def run_parsing_methods(self):
        """ runs all parsing methods to instantiate all attributes """

        # grabbing offset
        self.locate_offsets()

        # parsing header
        self.parse_sample_sheet_header()

        # parsing zumi
        self.parse_sample_sheet_zumi()

        # parsing diff exp
        self.parse_sample_sheet_diffexp()

        # parsing adapters
        self.parse_sample_sheet_adapters()

        # parsing cell data
        self.parse_sample_sheet_data()


    def create_adapter_whitelist(self):
        """ create barcode_whitelist text file """

        # need to direct toward specific path
        # will only exist in Docker
        self.cell_data.to_csv('barcode_white.txt', \
                                sep='\n', \
                                columns=['barcode_sequence'], \
                                header=False, \
                                index=False)


    def return_offsets(self):
        """ returns offset position info """
        return self.offset_pos

    def return_header_info(self):
        """ returns header info """
        return self.header_info

    def return_zumi_input(self):
        """ returns zumi input info """
        return self.zumi_input

    def return_diff_input(self):
        """ returns diff exp info """
        return self.diff_input

    def return_path_info(self):
        """ return path_info dict """
        return self.path_info

    def return_adapters(self):
        """ return adapter dict sequences """
        return self.adapter

    def return_cell_data(self):
        """ return cell barcode list """
        return self.cell_data



def main():
    """ run main for testing """
    pass

if __name__ == '__main__':
    main()
