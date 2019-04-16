#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

scrna_configure_pipeline.py

This is the main script that will take in the sample sheet and set the paths for other
various scripts and software to run throughout the pipeline. The script will
instantiate both the SampleSheetParser and ZumiConfigBuilder. The script will also
detect the amount of cores of the host server and designate cores as calculated.


Default run for STAR via zUMIs pipeline
STAR
    --genomeDir "STARidx"
    --runThreadN "p"
    --readFilesCommand zcat
    --sjdbGTFfile "gtf"
    --outFileNamePrefix "sample."
    --outSAMtype BAM Unsorted
    --outSAMmultNmax 1
    --outFilterMultimapNmax 50
    --outSAMunmapped Within
    --sjdbOverhang "readlength - 1"
    --twopassMode Basic
    --readFilesIn "cdnaread.filtered.fastq.gz"
additional_STAR_params: consider setting this as a default to handle many junctions
    --limitOutSJcollapsed 2000000
    --limitSjdbInsertNsj 2000000


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

CURRENT_DIR = os.getcwd()
LOG_DIR = 'logs'
DMEL_INDEX_DIR = 'dmel_star_idx_NOGTF'

if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)
    print("Made logs directory.")

if not os.path.exists(DMEL_INDEX_DIR):
    os.mkdir(DMEL_INDEX_DIR)
    print("Made star index directory.")


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



def samplesheet_zumi_build(sample_sheet_obj, zumi_config_obj):
    """ this creates the exchange of information from the SampleSheetParser
        object to the ZumiConfigBuilder object
    """

    # must run initially to set SampleSheetParser attributes
    sample_sheet_obj.run_parsing_methods()

    # creates the barcode whitelist text file for zumi
    sample_sheet_obj.create_adapter_whitelist(CURRENT_DIR)

    # returning all attribute information
    sam_header_info = sample_sheet_obj.return_header_info()
    sam_path_info = sample_sheet_obj.return_path_info()
    sam_zumi_info = sample_sheet_obj.return_zumi_input()
    sam_adapter_info = sample_sheet_obj.return_adapters()

    # need to make zumi output directory
    # ./basename_zumi_output
    zumi_config_obj.update_top_level(
        sam_header_info['run_name'] + '_' + sam_header_info['date'],
        CURRENT_DIR + '/' + sam_header_info['basename'] + '_zumi_output',
        16, # threads need to control for
        sam_zumi_info['zum_start_stage']
    )

    # passing trimmed fastqs paths to zumi config file
    zumi_config_obj.update_file_names(
        CURRENT_DIR + '/' + sam_path_info['trimmed_r1'],
        CURRENT_DIR + '/' + sam_path_info['trimmed_r2']
    )

    # need to better handle star index build path
    # setting reference paths for zumi output
    zumi_config_obj.update_reference_files(
        CURRENT_DIR + '/' + DMEL_INDEX_DIR,
        sam_path_info['annotation'] # gtf path
    )

    # passing filter cutoffs
    zumi_config_obj.update_filter_cutoffs(
        int(sam_zumi_info['bc_filter_num_bases']),
        int(sam_zumi_info['bc_filter_phred']),
        int(sam_zumi_info['umi_filter_num_bases']),
        int(sam_zumi_info['umi_filter_phred']),
    )

    # need to handle barcode whitelist path better
    zumi_config_obj.update_barcodes(
        CURRENT_DIR + '/barcode_whitelist.txt',
        int(sam_zumi_info['bc_ham_dist']),
    )

    # creates the nested dicts in the proper order to be dumped to YAML file
    # must run first before writing to yaml file
    zumi_config_obj.set_nested_dict()

    # contained in root directory of the repo
    zumi_config_yaml_path = zumi_config_obj.write_yaml(
        sam_header_info['basename'],
        CURRENT_DIR
    )

    # after all objects are parsed and yaml file is written, pass the path info
    # and adapter info back as dicts to properly run downstream scripts
    return sam_path_info, sam_adapter_info, zumi_config_yaml_path


def run_quality_control(**kwargs):
    """ takes kwargs and runs the qc wrapper function containing FASTQC and cutadapt

        params:
            threads : int for threading
            fastq_r1 : ab/path/to/raw_read_1.fastq, passed from SampleSheetParser
            fastq_r2 : ab/path/to/raw_read_2.fastq, passed from SampleSheetParser
            trimmed_r1 : ab/path/to/trimmed_read_1.fastq, must concat the current
                            working directory, or whichever desired direcotry,
                            string from SampleSheetParser only contains the filename
            trimmed_r2 : ab/path/to/trimmed_read_2.fastq, must concat the current
                            working directory, or whichever desired direcotry,
                            string from SampleSheetParser only contains the filename
            adapter_3 : str of 3' adapter sequence to trim, passed from SampleSheetParser
            adapter_5 : str of 5' adapter sequence to trim, passed from SampleSheetParser
    """

    # currently -m and -M are set to 12 for read1 and 20 for read 2
    # minimum length to trim to if a sequence must be trimmed
    run_qc_cmd = """python3 scrna_qc_pipeline.py
                    -t {threads}
                    -f {fastq_r1}
                    -p {fastq_r2}
                    -i {trimmed_r1}
                    -d {trimmed_r2}
                    -a {adapter_3}
                    -A {adapter_5}
                    -m 12
                    -M 20""".format(**kwargs)

    run_qc_formatted_args = shlex.split(run_qc_cmd)

    subprocess.run(
        run_qc_formatted_args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )


def build_star_index(**kwargs):
    """ takes kwargs and builds the STAR index

        params:
            threads : int for threading
            dmel_index_dir : output dir to put index files, needs to be created
                                prior to running STAR
            ref_fasta_file : absolute path to fasta file, grabs from SampleSheetParser
        returns:
            None
    """
    star_build_index_cmd = """STAR
                                --runMode genomeGenerate
                                --runThreadN {threads}
                                --genomeDir {dmel_index_dir}
                                --genomeFastaFiles {ref_fasta_file}
                            """.format(**kwargs)

    star_idx_formatted_args = shlex.split(star_build_index_cmd)

    subprocess.run(
        star_idx_formatted_args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )


def run_zumi_pipeline(zumi_yaml):
    """ runs the zumi pipeline. Steps include filtering, alignment, and counting

        params:
            zumi_yaml : path/to/zumi_config.yaml file constructed from the
                        ZumiConfigBuilder
        returns:
            None
    """

    zumi_cmd = """bash ./zUMIs-master.sh
                -y {zumi_config_yaml}""".format(zumi_config_yaml=zumi_yaml)

    zumi_formatted_args = shlex.split(zumi_cmd)

    subprocess.run(
        zumi_formatted_args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )


def main():
    """ parses sample sheet object, creates various config files, and tracks threads """

    args = arg_parser()

    # creates the representers needed for YAML
    setup_yaml()

    sample_sheet = args.sample_sheet
    LOGGER.info("Input args: %s", args)

    # intantiating each object
    sample_sheet_obj = SampleSheetParser(sample_sheet)
    zumi_config_obj = ZumiConfigBuilder()

    # parsing SampleSheetParser and exchanging information with ZumiConfigBuilder
    file_path_info, adapter_info, zumi_yaml = samplesheet_zumi_build(
        sample_sheet_obj,
        zumi_config_obj
    )

    # calls subprocess to run the scrna_qc_pipeline.py script
    run_quality_control(
        threads=16,
        fastq_read1=file_path_info['fastq_read1'],
        fastq_read2=file_path_info['fastq_read2'],
        trimmed_r1=CURRENT_DIR + file_path_info['trimmed_r1'],
        trimmed_r2=CURRENT_DIR + file_path_info['trimmed_r2'],
        adapter_3=adapter_info['adapter_3'],
        adapter_5=adapter_info['adapter_5']
    )

    # calls STAR by subprocess to build reference genome index
    build_star_index(
        threads=16,
        dmel_index_dir=CURRENT_DIR + '/' + DMEL_INDEX_DIR,
        ref_fasta_file=file_path_info['ref_genome']
    )

    # zumi yaml path is returned from samplesheet_zumi_build after writing the file
    run_zumi_pipeline(zumi_yaml)

if __name__ == '__main__':
    main()
