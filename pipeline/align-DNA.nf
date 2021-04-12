
nextflow.enable.dsl=2

log.info """\
   ===================================
   P I P E L I N E - A L I G N - D N A
   ===================================
   Boutros Lab

   Current Configuration:
   - input: 
      sample_name: ${params.sample_name}
      input_csv: ${params.input_csv}
      reference_fasta: ${params.reference_fasta}
      reference_fasta_index_files: ${params.reference_fasta_index_files}

   - output: 
      temp_dir: ${params.temp_dir}
      output_dir: ${params.output_dir}
      bam_output_dir: ${params.bam_output_dir}
      bam_output_filename: ${params.bam_output_filename}
      log_output_dir: ${params.log_output_dir}
      
   - options:
      save_intermediate_files = ${params.save_intermediate_files}
      cache_intermediate_pipeline_steps = ${params.cache_intermediate_pipeline_steps}
      max_number_of_parallel_jobs = ${params.max_number_of_parallel_jobs}
      blcds_registered_dataset_input = ${params.blcds_registered_dataset_input}
      blcds_registered_dataset_output = ${params.blcds_registered_dataset_output}

   Tools Used:
   - BWA-MEM2 and SAMtools: ${params.docker_image_bwa_and_samtools}
   - Picard Tools: ${params.docker_image_picardtools}
   - sha512sum: ${params.docker_image_sha512sum}
   - validate_params: ${params.docker_image_validate_params}

   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

include { validate_file } from './modules/validation.nf'
include { aligndna } from './modules/align-dna/workflow.nf'

workflow {
   Channel
      .fromPath(params.reference_fasta, checkIfExists: true)
      .set { ich_reference_fasta }

   Channel
      .fromPath(params.reference_fasta_index_files, checkIfExists: true)
      .set { ich_reference_index_files }

   // get the input fastq pairs
   Channel
      .fromPath(params.input_csv, checkIfExists: true)
      .splitCsv(header:true)
      .map { row -> 
         def read_group_name = "@RG" +
            "\\tID:" + row.read_group_identifier + ".Seq" + row.lane +
            "\\tCN:" + row.sequencing_center +
            "\\tLB:" + row.library_identifier +
            "\\tPL:" + row.platform_technology +
            "\\tPU:" + row.platform_unit +
            "\\tSM:" + row.sample

         // the library, sample and lane are used as keys downstream to group into 
         // sets of the same key for downstream merging
         return tuple(row.library_identifier,
            row.lane,
            read_group_name,
            row.read1_fastq,
            row.read2_fastq
            )
         }
      .set{ ich_samples }
   
   ich_samples
      .flatMap { library, lane, read_group_name, read1_fastq, read2_fastq ->
         [read1_fastq, read2_fastq]
         }
      .set { ich_samples_validate }

   validate_file(ich_samples_validate.mix(
      ich_reference_fasta,
      ich_reference_index_files
      ))

    aligndna(ich_samples,
       ich_reference_fasta,
       ich_reference_index_files
       )
   }
