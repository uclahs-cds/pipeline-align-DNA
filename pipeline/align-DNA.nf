
nextflow.enable.dsl=2

log.info """\
   ===================================
   P I P E L I N E - A L I G N - D N A
   ===================================
   Boutros Lab, version 7.0.4

   Current Configuration:
   - input: 
      sample_name: ${params.sample_name}
      input_csv: ${params.input_csv}
      reference_fasta_bwa: ${params.aligner.contains("BWA-MEM2") ? params.reference_fasta_bwa : "None"}
      reference_fasta_index_files_bwa: ${params.aligner.contains("BWA-MEM2") ? params.reference_fasta_index_files_bwa : "None"}
      reference_fasta_hisat2: ${params.aligner.contains("HISAT2") ? params.reference_fasta_hisat2 : "None"}
      reference_fasta_index_files_hisat2: ${params.aligner.contains("HISAT2") ? params.reference_fasta_index_files_hisat2 : "None"}

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
   - BWA-MEM2: ${params.aligner.contains("BWA-MEM2") ? params.docker_image_bwa_and_samtools : "None"}
   - HISAT2:  ${params.aligner.contains("HISAT2") ? params.docker_image_hisat2_and_samtools : "None"}
   - Picard Tools: ${params.docker_image_picardtools}
   - sha512sum: ${params.docker_image_sha512sum}
   - validate_params: ${params.docker_image_validate_params}

   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

include { validate_file } from './modules/run_validate.nf'
include { align_DNA_BWA_MEM2_workflow } from './modules/align_DNA_BWA_MEM2.nf'
include { align_DNA_HISAT2_workflow } from './modules/align_DNA_HISAT2.nf'

workflow {
   if (!(params.aligner.contains("BWA-MEM2") || params.aligner.contains("HISAT2"))) {
      throw new Exception('ERROR: Please specify at least one valid aligner! Options: BWA-MEM2, HISAT2')
      }

   // get the input fastq pairs
   Channel
      .fromPath(params.input_csv, checkIfExists: true)
      .splitCsv(header:true)
      .map { row -> 

         // the library, sample and lane are used as keys downstream to group into 
         // sets of the same key for downstream merging
         return tuple(row.library_identifier,
            row,
            row.lane,
            row.read1_fastq,
            row.read2_fastq
            )
         }
      .set{ ich_samples }
  
   ich_samples
      .flatMap { library, header, lane, read1_fastq, read2_fastq ->
         [read1_fastq, read2_fastq]
         }
      .set { ich_samples_validate }

   // Only create input channels for files which aligners are using
   if (params.aligner.contains("BWA-MEM2")) {
      Channel
         .fromPath(params.reference_fasta_bwa, checkIfExists: true)
         .set { ich_reference_fasta_bwa }
      Channel
         .fromPath(params.reference_fasta_index_files_bwa, checkIfExists: true)
         .set { ich_bwa_reference_index_files }
      align_DNA_BWA_MEM2_workflow(
         ich_samples,
         ich_samples_validate,
         ich_reference_fasta_bwa,
         ich_bwa_reference_index_files
         )
      }  
   if (params.aligner.contains("HISAT2")) {
      Channel
         .fromPath(params.reference_fasta_hisat2, checkIfExists: true)
         .set { ich_reference_fasta_hisat2 }
      Channel
         .fromPath(params.reference_fasta_index_files_hisat2, checkIfExists: true)
         .set { ich_hisat2_reference_index_files }
      align_DNA_HISAT2_workflow(
         ich_samples,
         ich_samples_validate,
         ich_reference_fasta_hisat2,
         ich_hisat2_reference_index_files
         )
      } 
   }
