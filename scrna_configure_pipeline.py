#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

scrna_configure_pipeline.py

This is the main script that will take in the sample sheet and set the paths for other
various scripts and software to run throughout the pipeline. The script will
instantiate both the SampleSheetParser and ZumiConfigBuilder. The script will also
detect the amount of cores of the host server and designate cores as calculated.
"""


import os
import shlex
import argparse
import subprocess
import logging
from collections import OrderedDict
import yaml
from yaml.resolver import Resolver
from samplesheet_parser import SampleSheetParser
from zumi_config_builder import ZumiConfigBuilder

LOG_DIR = 'logs'

if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)
    print("Made logs directory.")

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
FORMATTER = logging.Formatter('%(levelname)s:%(name)s:%(asctime)s:%(message)s')
FILE_HANDLER = logging.FileHandler("logs/main_config.log")
FILE_HANDLER.setFormatter(FORMATTER)
LOGGER.addHandler(FILE_HANDLER)


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Runs the single-cell RNA-seq pipeline.'
    )
    parser.add_argument('sample_sheet', type=str, help='Enter absolute/path/to/sample_sheet.csv')

    return parser.parse_args()


def represent_dictionary_order(self, dict_data):
    """ instantiates yaml dict mapping """
    return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())

def setup_yaml():
    """ adds the representer to the yaml instance """
    yaml.add_representer(OrderedDict, represent_dictionary_order)

    # remove resolver entries for On/Off/Yes/No
    for char_bool in "OoYyNn":
        if len(Resolver.yaml_implicit_resolvers[char_bool]) == 1:
            del Resolver.yaml_implicit_resolvers[char_bool]
        else:
            Resolver.yaml_implicit_resolvers[char_bool] = [
                x for x in Resolver.yaml_implicit_resolvers[char_bool]
                if x[0] != 'tag:yaml.org,2002:bool']


def build_qc_args(*args):
    """ builds arg dict for running qc wrapper; nonexhaustive list of args
        passable. will need to change this to include more if desired.

        args must be in in order of the keywords
    """

    # these args are not inclusive of all passable args
    # would need to fix this if other args are wanted
    qc_keywords = (
        'threads',
        'fastq_read1'
        'fastq_read2'
        'trimmed_r1'
        'trimmed_r2'
        'adapter_3'
        'adapter_5'
    )

    return dict(zip(qc_keywords, args))


def main():
    """ parses sample sheet object, creates various config files, and tracks threads """

    args = arg_parser()

    setup_yaml()

    LOGGER.info("Input args: %s", args)

    sample_sheet = args.sample_sheet

    LOGGER.info("Created SampleSheetParser Object")
    sample_sheet_obj = SampleSheetParser(sample_sheet)

    sample_sheet_obj.run_parsing_methods()
    LOGGER.info("Parsed sample sheet.")

    # creating barcode white list text file for zUMI
    sample_sheet_obj.create_adapter_whitelist()
    LOGGER.info("Created Barcode whitelist for zUMI.")

    LOGGER.info("Created ZumiConfigBuilder Object")
    zumi_config_obj = ZumiConfigBuilder()

    # dict contains all relevant file paths
    file_path_info = sample_sheet_obj.return_path_info()

    # dict contains adapter trimming sequences
    adapter_info = sample_sheet_obj.return_adapters()

    qc_kwarg = build_qc_args(
        16,
        file_path_info['fastq_read1'],
        file_path_info['fastq_read2'],
        file_path_info['trimmed_r1'],
        file_path_info['trimmed_r2'],
        adapter_info['adapter_3'],
        adapter_info['adapter_5']
    )

    # need to add kwarg dict to initiate run
    run_qc_cmd = """python3 scrna_qc_pipeline.py
                    -t {threads}
                    -f {fastq_r1}
                    -p {fastq_r2}
                    -i {trimmed_r1}
                    -d {trimmed_r2}
                    -a {adapter_3}
                    -A {adapter_5}
                    -m 12
                    -M 20""".format(**qc_kwarg)

    run_qc_formatted_args = shlex.split(run_qc_cmd)

    subprocess.run(run_qc_formatted_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

if __name__ == '__main__':
    main()
