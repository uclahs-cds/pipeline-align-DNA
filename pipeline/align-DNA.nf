
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
      reference_fasta_index_files: ${params.aligner.contains("BWA-MEM2") ? params.bwa_reference_fasta_index_files : ""} ${params.aligner.contains("HISAT2") ? params.hisat2_reference_fasta_index_files : ""}

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
   - ${params.aligner}: ${params.aligner.contains("BWA-MEM2") ? params.docker_image_bwa_and_samtools : ""} ${params.aligner.contains("HISAT2") ? params.docker_image_hisat2_and_samtools : ""}
   - Picard Tools: ${params.docker_image_picardtools}
   - sha512sum: ${params.docker_image_sha512sum}
   - validate_params: ${params.docker_image_validate_params}

   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

include { validate_file; validate_file as validate_bwa_output_file; validate_file as validate_hisat2_output_file } from './modules/run_validate.nf'
include { align_DNA_BWA_MEM2 } from './modules/align_DNA_BWA_MEM2.nf'
include { align_DNA_HISAT2 } from './modules/align_DNA_HISAT2.nf'
include { PicardTools_SortSam; PicardTools_SortSam as PicardTools_SortSam_HISAT2 } from './modules/sort_bam_picardtools.nf'
include { PicardTools_MarkDuplicates; PicardTools_MarkDuplicates as PicardTools_MarkDuplicates_HISAT2 } from './modules/mark_duplicate_picardtools.nf'
include { PicardTools_BuildBamIndex; PicardTools_BuildBamIndex as PicardTools_BuildBamIndex_HISAT2 } from './modules/index_bam_picardtools.nf'
include { Generate_Sha512sum; Generate_Sha512sum as Generate_Sha512sum_HISAT2 } from './modules/check_512sum.nf'

workflow {
   Channel
      .fromPath(params.reference_fasta, checkIfExists: true)
      .set { ich_reference_fasta }
   
   Channel
      .fromPath(params.bwa_reference_fasta_index_files, checkIfExists: params.aligner.contains("BWA-MEM2"))
      .set { ich_bwa_reference_index_files }

   Channel
      .fromPath(params.hisat2_full_reference_fasta_index_files, checkIfExists: params.aligner.contains("HISAT2"))
      .set { ich_hisat2_reference_index_files }

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

   if (params.aligner.contains("BWA-MEM2") && params.aligner.contains("HISAT2")) {
      validate_file(ich_samples_validate.mix(
         ich_reference_fasta,
         ich_bwa_reference_index_files,
         ich_hisat2_reference_index_files
         ))
      }
   else if (params.aligner.contains("BWA-MEM2")) {
      validate_file(ich_samples_validate.mix(
         ich_reference_fasta,
         ich_bwa_reference_index_files
         ))
      }
   else if (params.aligner.contains("HISAT2")) {
      validate_file(ich_samples_validate.mix(
         ich_reference_fasta,
         ich_hisat2_reference_index_files
         ))
      }

   if (params.aligner.contains("BWA-MEM2")) {
      aligner_output_dir = "${params.bam_output_dir}-BWA-MEM2"
      align_DNA_BWA_MEM2(
         ich_samples,
         ich_reference_fasta,
         ich_bwa_reference_index_files.collect()
         )
      PicardTools_SortSam(align_DNA_BWA_MEM2.out.bam)
      PicardTools_MarkDuplicates(PicardTools_SortSam.out.bam.collect(), aligner_output_dir)
      PicardTools_BuildBamIndex(PicardTools_MarkDuplicates.out.bam, aligner_output_dir)
      Generate_Sha512sum(PicardTools_BuildBamIndex.out.bai.mix(PicardTools_MarkDuplicates.out.bam), aligner_output_dir)
         validate_bwa_output_file(
            PicardTools_MarkDuplicates.out.bam.mix(
               PicardTools_BuildBamIndex.out.bai,
               Channel.from(params.temp_dir, params.output_dir)
               )
            )
      }  
   if (params.aligner.contains("HISAT2")) {
      aligner_output_dir = "${params.bam_output_dir}-HISAT2" 
      align_DNA_HISAT2(
         ich_samples,
         ich_reference_fasta,
         ich_hisat2_reference_index_files.collect()
         )
      PicardTools_SortSam_HISAT2(align_DNA_HISAT2.out.bam)   
      PicardTools_MarkDuplicates_HISAT2(PicardTools_SortSam_HISAT2.out.bam.collect(), aligner_output_dir)
      PicardTools_BuildBamIndex_HISAT2(PicardTools_MarkDuplicates_HISAT2.out.bam, aligner_output_dir)
      Generate_Sha512sum_HISAT2(PicardTools_BuildBamIndex_HISAT2.out.bai.mix(PicardTools_MarkDuplicates_HISAT2.out.bam), aligner_output_dir)
         validate_hisat2_output_file(
            PicardTools_MarkDuplicates_HISAT2.out.bam.mix(
               PicardTools_BuildBamIndex_HISAT2.out.bai,
               Channel.from(params.temp_dir, params.output_dir)
               )
            )
      } 
   }
