"""write_dna_align_config_file
Usage:
    write_dna_align_config_file.py param
    write_dna_align_config_file.py example
    write_dna_align_config_file.py <input_file> (bwa-mem2 | hisat2) <reference_fasta>
        <output_dir> <temp_dir>
        [--save_intermediate_files]
        [--cache_intermediate_pipeline_steps]
        [--HISAT2_fasta_index=<path>]
        [--blcds_registered_dataset_input]
        [--blcds_registered_dataset_output]
        [--save_bam_and_log_to_blcds]
        [--blcds_disease_id=<disease_id>]
        [--tblcds_dataset_id=<dataset_id>]
        [--tblcds_patient_id=<patient_id>]
        [--tblcds_sample_id=<sample_id>]
        [--tblcds_analyte=<analyte>]
        [--tblcds_technology=<technology>]
        [--tblcds_mount_dir=</data>]

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
        (bwa-mem2 | hisat2): choose between 2 aligners (for HISAT2 fasta_index option is required)
        <reference_fasta>: absolute path to the refernce fasta file
        <output_dir>: absolute path to the output directory
        <temp_dir>: absolute path to a temp directory

    Options:
        --HISAT2_fasta_index=<path>: Only if hisat2 is chosen, the reference fasta index should be given

        --save_intermediate_files: use to save intermediate files. [default: None]
        --cache_intermediate_pipeline_steps default  : [default: None]
        --blcds_registered_dataset_input: use if the data input fastq files are registered 
            in the Boutros Lab [default: None]
        --blcds_registered_dataset_output  : use to redirect output files directly to the Boutros 
            Lab data storage [default: None]

        --save_bam_and_log_to_blcds: use in order to save output bam and log directly to blcds data 
            storage [default: None]
        if save_bam_and_log_to_blcds is used the following parameters must be given
            disease_id (by using --blcds_disease_id)
            dataset_id (by using --tblcds_dataset_id)
            patient_id (by using --tblcds_patient_id)
            sample_id (by using --tblcds_sample_id)
            analyte (for example: DNA or RNA) (by using --tblcds_analyte)
            technology (for example: WGS or WTS (by using --tblcds_technology)
            mount_dir (by using --tblcds_mount_dir)
    """
    print(params)


def print_example():
    """
    printing commands for example
    """
    example = "example command 1: here the bwa-mem2 aligner is chosen:\n\n" \
        + "python3 write_DNA_align_config_file.py" \
        + "path/to/input/name_of_sample_with_no_spaces bwa-mem2 /path/to/bwa/fasta/genome.fa" \
        + "where/to/save/outputs/ /local/disk/for/temp/file/dir/\n\n" \
        + "example command 2: here the hisat2 aligner is chosen and " \
        + "--cache_intermediate_pipeline_steps, and --save_bam_and_log_to_blcds flages are used:" \
        + "\n\npython3 write_dna_align_config_file.py" \
        + "path/to/input/name_of_sample_with_no_spaces hisat2 /path/to/bwa/fasta/genome.fa" \
        + "where/to/save/outputs/ /local/disk/for/temp/file/dir/ --save_bam_and_log_to_blcds" \
        + "--cache_intermediate_pipeline_steps --HISAT2_fasta_index=path/to/index" \
        + "--blcds_disease_id=my_disease_id --tblcds_patient_id=my_patient_id" \
        + "--tblcds_sample_id=my_sample_id --tblcds_analyte=RNA --tblcds_technology=WTS" \
        + "--tblcds_mount_dir=/hot --tblcds_dataset_id=my_dataset-i"
    print(example)

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
            'hisat2_index_prefix = "' + args['--HISAT2_fasta_index'] + '"\n\t' \
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
        # test if all required flages were used
        if not(args['--blcds_disease_id']
            and args['--tblcds_dataset_id']
            and args['--tblcds_patient_id']
            and  args['--tblcds_sample_id']
            and args['--tblcds_analyte']
            and args['--tblcds_technology']
            and args['--tblcds_mount_dir']):
            raise TypeError(
                "if save_bam_and_log_to_blcds is used, disease_id, dataset_id, patient_id" \
                "sample_id, analyte, technology', and mount_dir must be given"
                )
        return '\tblcds_cluster_slurm = true\n' \
            + '\tblcds_disease_id = "' + args['--blcds_disease_id'] + '"\n' \
            + '\tblcds_dataset_id = "' + args['--tblcds_dataset_id'] + '"\n' \
            + '\tblcds_patient_id = "' + args['--tblcds_patient_id'] + '"\n' \
            + '\tblcds_sample_id  = "' + args['--tblcds_sample_id'] + '"\n' \
            + '\tblcds_analyte = "' + args['--tblcds_analyte'] + '"\n' \
            + '\tblcds_technology = "' + args['--tblcds_technology'] + '"\n' \
            + '\tblcds_mount_dir = "' + args['--tblcds_mount_dir'] + '"\n'
    return ''


if __name__ == '__main__':
    arguments = docopt(__doc__)

    if arguments['param']:
        print_params()
        sys.exit()

    if arguments['example']:
        print_example()
        sys.exit()

    sample_name = os.path.basename(arguments['<input_file>'])
    # create the string for writing to the config file
    string_for_writing = "final String SAMPLE = '" + sample_name + "'\n" \
        + 'params {' + '\n' \
        + '\t' + 'sample_name = SAMPLE' + '\n' \
        + '\t' + 'input_csv = "' \
        + os.path.join(os.path.dirname(arguments['<input_file>']), '${SAMPLE}.csv') + '"\n' \
        + '\t' + aligner_type_choice(arguments) + '\n' \
        + '\t' + 'output_dir = "' + os.path.join(arguments['<output_dir>'], '${SAMPLE}') + '"\n' \
        + '\t' + 'temp_dir = "' + arguments['<temp_dir>'] + '"\n' \
        + '\t' + is_save_inter(arguments) + '\n' \
        + '\t' + is_cache_intermediate_pipeline_steps(arguments) + '\n' \
        + '\t' + is_blcds_registered_dataset_input(arguments) + '\n' \
        + '\t' + blcds_registered_dataset_output(arguments) + '\n' \
        + is_save_bam_and_log_to_blcds(arguments) \
        + '\t' + '}' + '\n' \
        + 'includeConfig "${projectDir}/config/methods.config"\n' \
        + 'methods.setup()'
    # writing to the file
    with open(sample_name + '_DNA_align.config', 'w') as config:
        config.write(string_for_writing)
