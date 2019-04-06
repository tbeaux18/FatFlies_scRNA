#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

zumi_config_builder.py

Must include yaml representer to preserve order
zUMI requires all paths to be absolute; will need to take this into account
for creating the barcode text file and trimmed fastq samples
"""


from collections import OrderedDict
import yaml

def represent_dictionary_order(self, dict_data):
    """ instantiates yaml dict mapping """
    return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())

def setup_yaml():
    """ adds the representer to the yaml instance """
    yaml.add_representer(OrderedDict, represent_dictionary_order)


class ZumiConfigBuilder:
    """
    ZumiConfigBuilder object instantiates the appropriate dictionaries for
    easy value updates in the format provided by the zUMIs pipeline.

    Order matters for how the yaml file is produced. This is a hacked way of
    getting the proper format for the correct structure so zUMIs will work.

    Attributes:
        top_yaml_dict
        sequence_file_dict
        reference_dict
        filter_cutoff_dict
        barcode_dict
        counting_opts_dict

    Instance Methods:
        set_nested_dict
        return_top_yaml
        return_sequence_dict
        return_reference_dict
        return_filter_dict
        return_barcode_dict
        return_count_dict
        write_2_yaml

    """

    def __init__(self):

        self.top_yaml_dict = OrderedDict([
            ('project', None),
            ('sequence_files', None),
            ('reference', None),
            ('out_dir', None),
            ('num_threads', None),
            ('mem_limit', None),
            ('filter_cutoffs', None),
            ('barcodes', None),
            ('counting_opts', None),
            ('make_stats', 'yes',),
            ('which_Stage', None),
            ('samtools_exec', 'samtools'),
            ('Rscript_exec', 'Rscript'),
            ('STAR_exec', 'STAR'),
            ('pigz_exec', 'pigz')
        ])

        self.sequence_file_dict = OrderedDict([
            ('file1', OrderedDict([
                ('name', None), # update from SampleSheetParser
                ('base_definition', ['BC(7-12)', 'UMI(1-6)'])
            ])
            ),
            ('file2', OrderedDict([
                ('name', None), # update from SampleSheetParser
                ('base_definition', ['cDNA(1-50)'])
            ])
            )
        ])

        self.reference_dict = OrderedDict([
            ('STAR_index', None), # update from SampleSheetParser
            ('GTF_file', None), # update from SampleSheetParser
            ('additional_files', None),
            ('additional_STAR_params', None)
        ])

        self.filter_cutoff_dict = OrderedDict([
            ('BC_filter', OrderedDict([
                ('num_bases', 1), # update from SampleSheetParser
                ('phred', 20) # update from SampleSheetParser
            ])
            ),
            ('UMI_filter', OrderedDict([
                ('num_bases', 1), # update from SampleSheetParser
                ('phred', 20) # update from SampleSheetParser
            ])
            )
        ])

        self.barcode_dict = OrderedDict([
            ('barcode_num', None),
            ('barcode_file', None), # file path from SampleSheetParser
            ('automatic', 'yes'),
            ('BarcodeBinning', None), # int from SampleSheetParser
            ('nReadsperCell', 100)
        ])

        self.counting_opts_dict = OrderedDict([
            ('introns', 'yes'),
            ('downsampling', 0),
            ('strand', 0),
            ('Ham_Dist', 1), # need to include this as an option on sample_sheet
            ('velocyto', 'no'),
            ('primaryHit', 'yes'),
            ('twoPass', 'yes')
        ])

    def set_nested_dict(self):
        """ setting nested keys with their respective values """

        # these set the nested dicts needed for proper yaml structure
        self.top_yaml_dict['sequence_files'] = self.sequence_file_dict
        self.top_yaml_dict['reference'] = self.reference_dict
        self.top_yaml_dict['filter_cutoffs'] = self.filter_cutoff_dict
        self.top_yaml_dict['barcodes'] = self.barcode_dict
        self.top_yaml_dict['counting_opts'] = self.counting_opts_dict

    def return_top_yaml(self):
        """ return the main yaml dict """
        return self.top_yaml_dict

    def return_sequence_dict(self):
        """ returns sequence file dict """
        return self.sequence_file_dict

    def return_reference_dict(self):
        """ returns reference path dicts """
        return self.reference_dict

    def return_filter_dict(self):
        """ returns filter cutoff dict """
        return self.filter_cutoff_dict

    def return_barcode_dict(self):
        """ returns barcode opts dict """
        return self.barcode_dict

    def return_count_dict(self):
        """ returns counting opts dict """
        return self.counting_opts_dict

    def write_2_yaml(self):
        """ writes the top level yaml dict to yaml file """
        with open('test-zumi.yaml', 'w') as outfile:
            yaml.dump(self.top_yaml_dict, outfile, default_flow_style=False, default_style=None)



def main():
    """ run main for testing """

    setup_yaml()

    z = ZumiConfigBuilder()

    z.set_nested_dict()

    z.write_2_yaml()

if __name__ == '__main__':
    main()
