"""write_dna_align_config_file
Usage:
    write_dna_align_config_file.py param
    write_dna_align_config_file.py <input_file> (bwa-mem2 | hisat2) <reference_fasta>
        <output_dir> <temp_dir> [--save_intermediate_files] [--cache_intermediate_pipeline_steps]
        [--HISAT2_fasta_index=<path>] [--blcds_registered_dataset_input]
        [--blcds_registered_dataset_output] [--save_bam_and_log_to_blcds]

Options:
  -h --help           Show this screen.
  run 'write_dna_align_config_file.py param' to learn about options
"""

### this script creates a config file that is required to run the align DNA pipeline ###############

import sys
import os
from docopt import docopt

def print_params():
    """
    printing info about the paramters and options
    """
    params = """
    Parameters:
        <input_file>: absolute path to the input fastq file
        (BWA-MEM2 | HISAT2): choose between 2 aligners (for HISAT2 fasta_index option is required)
        <reference_fasta>: absolute path to the refernce fasta file
        <output_dir>: absolute path to the output directory
        <temp_dir>: absolute path to a temp directory
    Options:
        --HISAT2_fasta_index=<path>: Only if HISAT2 is chosen, the reference fasta index should be given
        --save_intermediate_files: enter to save intermediate files. [default: None]
        --cache_intermediate_pipeline_steps default  : [default: None]
        --blcds_registered_dataset_input: set to true if the data input fastq files are registered 
            in the Boutros Lab [default: None]
        --blcds_registered_dataset_output  : set to true to redirect output files directly to the Boutros 
            Lab data storage [default: None]
        --save_bam_and_log_to_blcds: enter in order to save output bam and log directly to blcds data 
            storage [default: None]
    """
    print(params)

def aligner_type_choice(args):
    """
    :param args: docopt arguments
    :return: the lines to add the file depending on the type of chosen aligner
    """
    if args['bwa-mem2']:
        return  'reference_fasta_bwa = "' + args['<reference_fasta>'] + '"\n\t' \
            + 'aligner = ["BWA-MEM2"]'
    if args['hisat2']:
        if args['--HISAT2_fasta_index'] is None:
            raise TypeError(
                "for HISAT2 aligner, --HISAT2_fasta_index=<path> must be given' "
                )
        return 'reference_fasta_hisat2 = "' + args['<reference_fasta>'] + '"\n\t' +  \
            'reference_fasta_index_prefix = "' + args['--HISAT2_fasta_index'] + '"\n\t' \
            + 'aligner = ["HISAT2"]'
    return ''
def is_save_inter(args):
    """
    :param args: docopt arguments
    :return: if save_intermediate_files flag is used, return 'the line that should be printed =
    true', otherwise false
    """
    if args['--save_intermediate_files']:
        return 'save_intermediate_files = true'
    return 'save_intermediate_files = false'

def is_cache_intermediate_pipeline_steps(args):
    """
    :param args: docopt arguments
    :return: if cache_intermediate_pipeline_steps flag is used, return 'the line that should be
    printed = true', otherwise false
    """
    if args['--cache_intermediate_pipeline_steps']:
        return 'cache_intermediate_pipeline_steps = true'
    return 'cache_intermediate_pipeline_steps = false'

def is_blcds_registered_dataset_input(args):
    """
    :param args: docopt arguments
    :return: if blcds_registered_dataset_input flag is used, return 'the line that should be
    printed = true', otherwise false
    """
    if args['--blcds_registered_dataset_input']:
        return 'blcds_registered_dataset_input = true'
    return 'blcds_registered_dataset_input = false'

def blcds_registered_dataset_output(args):
    """
    :param args: docopt arguments
    :return: if blcds_registered_dataset_output flag is used, return 'the line that should be
    printed = true', otherwise false
    """
    if args['--blcds_registered_dataset_output']:
        return 'blcds_registered_dataset_output = true'
    return 'blcds_registered_dataset_output = false'
def is_save_bam_and_log_to_blcds(args):
    """
    :param args: docopt arguments
    :return: if save_bam_and_log_to_blcds flag is used, return several lines that should be
        writing to the file. If not, nothing is printed (return '')
    """
    if args['--save_bam_and_log_to_blcds']:
        return '\tblcds_cluster_slurm = true\n' \
            + '\tblcds_disease_id = "disease_id"\n' \
            + '\tblcds_dataset_id = "dataset_id"\n' \
            + '\tblcds_patient_id = "patient_id"\n' \
            + '\tblcds_sample_id  = "sample_id"\n' \
            + '\tblcds_analyte = "DNA"\n' \
            + '\tblcds_technology = "WGS"\n' \
            + '\tblcds_mount_dir = "/data"\n'
    return ''


if __name__ == '__main__':
    arguments = docopt(__doc__)

    if arguments['param']:
        print_params()
        sys.exit()

    sample_name = os.path.basename(arguments['<input_file>'])

    # create the string for writing to the config file
    string_for_writing = "final String SAMPLE = '" + sample_name + "'\n" \
        + 'params {' + '\n' \
        + '\t' + 'sample_name = SAMPLE' + '\n' \
        + '\t' + 'input_csv = "' + os.path.dirname(arguments['<input_file>']) + '${SAMPLE}.csv' \
        + '"\n' \
        + '\t' + aligner_type_choice(arguments) + '\n' \
        + '\t' + 'output_dir = "' + os.path.join(arguments['<output_dir>'], '${SAMPLE}') + '"\n' \
        + '\t' + 'temp_dir = "' + arguments['<temp_dir>'] + '"\n' \
        + '\t' + is_save_inter(arguments) + '\n' \
        + '\t' + is_cache_intermediate_pipeline_steps(arguments) + '\n' \
        + '\t' + is_blcds_registered_dataset_input(arguments) + '\n' \
        + '\t' + blcds_registered_dataset_output(arguments) + '\n' \
        + is_save_bam_and_log_to_blcds(arguments) \
        + 'includeConfig "${projectDir}/config/methods.config\n' \
        + 'methods.setup()'
    # writing to the file
    with open('DNA_align_config_file', 'w') as config:
        config.write(string_for_writing)
