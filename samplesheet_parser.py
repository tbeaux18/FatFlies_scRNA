#!/usr/bin/env python3
"""
@author: Timothy Baker

samplesheet_parser.py

"""

import os
import pandas as pd


class SampleSheetParser:
    """

    A SampleSheetParser object that parses the input sample sheet in specific
    format and holds all instantiated paths.

    params:
        sample_sheet : ab/path/to/sample_sheet.csv

    attributes:
        fastq_read1
        fastq_read2
        ref_genome
        annotation
        adapter_3
        adapter_5
        barcode_seq


    """

    def __init__(self, sample_sheet):
        self.sample_sheet = sample_sheet
        self.fastq_read1 = None
        self.fastq_read2 = None
        self.ref_genome = None
        self.annotation = None
        self.adapter_3 = None
        self.adapter_5 = None
        self.cell_data = None


    def parse_sample_sheet(self):
        """ parsing sample sheet and instantiating attributes """

        with open(self.sample_sheet, 'r') as csv_handle:

            line = csv_handle.readline()

            while line:

                if line.startswith('[HEADER]'):
                    header_offset = csv_handle.tell()
                if line.startswith('[SETTINGS]'):
                    settings_offset = csv_handle.tell()
                if line.startswith('[DATA]'):
                    data_offset = csv_handle.tell()

                line = csv_handle.readline()

            # create header byte load
            byte_load = (settings_offset - header_offset) - 28

            # change position to header offset
            csv_handle.seek(header_offset)
            for line in csv_handle.readlines(byte_load):
                line_lst = line.split(',')
                if line_lst[0] == 'Fastq_Read1':
                    self.fastq_read1 = line_lst[1]
                if line_lst[0] == 'Fastq_Read2':
                    self.fastq_read2 = line_lst[1]
                if line_lst[0] == 'Ref_genome':
                    self.ref_genome = line_lst[1]
                if line_lst[0] == 'Annotation':
                    self.annotation = line_lst[1]

            # changing file position to adapter offset
            csv_handle.seek(settings_offset)

            # setting adapter attributes
            adapter_list = csv_handle.readline().split(',')
            self.adapter_3 = adapter_list[1]
            self.adapter_5 = adapter_list[2]

            # changing file position to data offset to load into dataframe
            csv_handle.seek(data_offset)

            # pandas dataframe for easy
            self.cell_data = pd.read_csv(csv_handle)


    def create_adapter_whitelist(self):
        """ create barcode_whitelist text file """
        self.cell_data.to_csv('barcode_white.txt', \
                                sep='\n', \
                                columns=['barcode_sequence'], \
                                header=False, \
                                index=False)
        
    def create_design_file(self):
        """ create design_matrix text file """
        self.cell_data.to_csv('design_matrix.txt', \
                                sep="\n", \
                                columns=['sample_name','barcode_sequence','experiment_group'],
                                header=False,
                                index=False)
                          

    def return_read1(self):
        """ return fastq_read1 path """
        return self.fastq_read1

    def return_read2(self):
        """ return fastq_read2 path """
        return self.fastq_read2

    def return_genome(self):
        """ return reference genome path """
        return self.ref_genome

    def return_annotation(self):
        """ return annotation gtf file path """
        return self.annotation

    def return_adapter_3(self):
        """ return 3' adapter sequence """
        return self.adapter_3

    def return_adapter_5(self):
        """ return 5' adapter sequence """
        return self.adapter_5

    def return_barcode_seq(self):
        """ return cell barcode list """
        return self.cell_data



def main():
    """ run main for testing """

    sample_object = SampleSheetParser('scrna_pipeline_samplesheet_template.csv')

    sample_object.parse_sample_sheet()

    sample_object.create_adapter_whitelist()

    sample_object.create_design_file()

if __name__ == '__main__':
    main()
