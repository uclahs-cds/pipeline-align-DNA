#!/usr/bin/env nextflow 

def dockeri_BWA_and_SAMtools = "blcdsdockerregistry/align-dna:${params.bwa_version}_samtools-1.10"
def dockeri_validate_params = "blcdsdockerregistry/validate:1.0.0"

// output details of the pipeline run to stdout
log.info """\
   ===========================================================
   P I P E L I N E - G E N E R A T E - G E N O M E - I N D E X
   ===========================================================
   Boutros Lab

   Current Configuration:
   - input:
        reference_fasta:    ${params.reference_fasta}

   - output:
        output_dir:         ${params.output_dir}

   Tools Used:
   - BWA-MEM2 and SAMtools: ${dockeri_BWA_and_SAMtools}

   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

Channel.fromPath(params.reference_fasta).into { ich_reference_fasta; ich_reference_fasta_validate }

process validate_inputs {
    container dockeri_validate_params

    input:
        path(file_to_validate) from ich_reference_fasta_validate

    script:
    """
    set -euo pipefail

    python -m validate -t file-input ${file_to_validate}
    """
}

process generate_index {
    container dockeri_BWA_and_SAMtools

    publishDir path: params.output_dir,
        pattern: "${reference_fasta.getName()}.*",
        mode: 'copy'

    publishDir path: "${params.log_output_dir}/generate_index",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "log${file(it).getName()}" }

    // use "each" so the the reference files are passed through for each fastq pair alignment 
    input: 
        file(reference_fasta) from ich_reference_fasta

    // output the lane information in the file name to differentiate bewteen aligments of the same
    // sample but different lanes
    output:
        file "${reference_fasta.getName()}.*" into och_generate_index
        file ".command.*"

    script:
    """
    set -euo pipefail
    bwa-mem2 index ${reference_fasta}
    """
}

process validate_outputs {
   container dockeri_validate_params

   input:
      path(file_to_validate) from och_generate_index.flatten()

   script:
   """
   set -euo pipefail
   python -m validate -t file-input ${file_to_validate}
   """
}
