
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
   - Aligner: ${params.docker_image_selected}
   - Picard Tools: ${params.docker_image_picardtools}
   - sha512sum: ${params.docker_image_sha512sum}
   - validate_params: ${params.docker_image_validate_params}

   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

include { validate_file; validate_file as validate_output_file } from './modules/run_validate.nf'
include { align_DNA_BWA_MEM2 } from './modules/align_DNA_BWA_MEM2.nf'
include { align_DNA_HISAT2 } from './modules/align_DNA_HISAT2.nf'
include { PicardTools_SortSam } from './modules/sort_bam_picardtools.nf'
include { PicardTools_MarkDuplicates } from './modules/mark_duplicate_picardtools.nf'
include { PicardTools_BuildBamIndex } from './modules/index_bam_picardtools.nf'
include { Generate_Sha512sum } from './modules/check_512sum.nf'

workflow {
   Channel
      .fromPath(params.reference_fasta, checkIfExists: true)
      .set { ich_reference_fasta }

   ref_fasta_index_files = params.reference_fasta_index_files
   if (params.aligner == "HISAT2") {
      ref_fasta_index_files = ref_fasta_index_files + ".*.ht2"
      }
   
   Channel
      .fromPath(ref_fasta_index_files, checkIfExists: true)
      .set { ich_reference_index_files }

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

   validate_file(ich_samples_validate.mix(
      ich_reference_fasta,
      ich_reference_index_files
      ))

   if (params.aligner == "BWA-MEM2") {
      align_DNA_BWA_MEM2(
         ich_samples,
         ich_reference_fasta,
         ich_reference_index_files.collect()
         )
      PicardTools_SortSam(align_DNA_BWA_MEM2.out.bam)
      }  
   else if (params.aligner == "HISAT2") { 
      align_DNA_HISAT2(
         ich_samples,
         ich_reference_fasta,
         ich_reference_index_files.collect()
         )
      PicardTools_SortSam(align_DNA_HISAT2.out.bam)   
      } 
   PicardTools_MarkDuplicates(PicardTools_SortSam.out.bam.collect())
   PicardTools_BuildBamIndex(PicardTools_MarkDuplicates.out.bam)
   Generate_Sha512sum(PicardTools_BuildBamIndex.out.bai.mix(PicardTools_MarkDuplicates.out.bam))
   validate_output_file(
      PicardTools_MarkDuplicates.out.bam.mix(
         PicardTools_BuildBamIndex.out.bai,
         Channel.from(params.temp_dir, params.output_dir)
         )
      )

   }
