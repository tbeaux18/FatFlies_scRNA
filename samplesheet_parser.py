#!/usr/bin/env python3
"""
@author: Timothy Baker

samplesheet_parser.py

"""

import sys
import csv


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
        self.barcode_seq = []
        # create dict attribute that holds info after data
        
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
        return self.barcode_seq








def main():
    pass

if __name__ == '__main__':
    main()
