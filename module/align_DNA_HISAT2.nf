// The follwing process runs both alignment and SAM conversion in the same process.
// While this is normally considered to go against the best practices for processes,
// here it actually saves cost, time, and memory to directly pipe the output into
// samtools due to the large size of the uncompressed SAM files.
include { generate_standard_filename } from '../external/nextflow-modules/modules/common/generate_standardized_filename/main.nf'
include { run_sort_SAMtools ; run_merge_SAMtools} from './samtools.nf'
include {
   run_validate_PipeVal as validate_input_HISAT2
   run_validate_PipeVal as validate_output_HISAT2
   } from '../external/nextflow-modules/modules/PipeVal/validate/main.nf' addParams(
         options: [
            log_output_dir: "${params.log_output_dir}/process-log/${params.hisat2_version}",
            docker_image_version: params.pipeval_version,
            main_process: "./"
            ]
      )
include { run_MarkDuplicate_Picard } from './mark_duplicate_picardtools.nf'
include { run_MarkDuplicatesSpark_GATK } from './mark_duplicates_spark.nf'
include { generate_sha512sum } from './check_512sum.nf'
include { remove_intermediate_files } from '../external/nextflow-modules/modules/common/intermediate_file_removal/main.nf' addParams(
   options: [
      save_intermediate_files: params.save_intermediate_files,
      output_dir: params.output_dir_base,
      log_output_dir: "${params.log_output_dir}/process-log/${params.hisat2_version}"
      ]
   )

process align_DNA_HISAT2 {
   container params.docker_image_hisat2_and_samtools
   publishDir path: "${params.output_dir_base}/${params.hisat2_version}/intermediate/${task.process.split(':')[1].replace('_', '-')}",
      enabled: params.save_intermediate_files,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: "${params.log_output_dir}/process-log/${params.hisat2_version}/${task.process.split(':')[1].replace('_', '-')}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${library}/${lane}/log${file(it).getName()}" }

   // use "each" so the the reference files are passed through for each fastq pair alignment
   input:
      tuple(val(library),
         val(header),
         val(lane),
         path(read1_fastq),
         path(read2_fastq)
         )
      each path(ref_fasta)
      path(ich_reference_index_files)

   // output the lane information in the file name to differentiate bewteen aligments of the same
   // sample but different lanes
   output:
      tuple val(library),
         val(lane),
         path("${lane_level_bam}"), emit: bam
      path(".command.*")

   script:

   lane_level_bam = generate_standard_filename(params.hisat2_version, params.dataset_id, params.sample_id, [additional_information: "${library}-${lane}.bam"])

   """
   set -euo pipefail

   hisat2 \
      -p ${task.cpus} \
      --rg-id ${header.read_group_identifier}.Seq${header.lane} \
      --rg LB:${header.library_identifier} \
      --rg PL:${header.platform_technology} \
      --rg PU:${header.platform_unit} \
      --rg SM:${header.sample} \
      --no-spliced-alignment \
      -x ${file(ich_reference_index_files.get(0).baseName).baseName} \
      -1 ${read1_fastq} \
      -2 ${read2_fastq} | \
   samtools \
      view \
      -@ ${task.cpus} \
      -S \
      -b > \
      ${lane_level_bam}
   """
   }

workflow align_DNA_HISAT2_workflow {
   aligner_output_dir = (params.ucla_cds_registered_dataset_output) ? "${params.output_dir_base}/${params.hisat2_version}/BAM-${params.hisat2_uuid}" : "${params.output_dir_base}/${params.hisat2_version}/output"
   aligner_intermediate_dir = (params.ucla_cds_registered_dataset_output) ? "${params.output_dir_base}/${params.hisat2_version}/BAM-${params.hisat2_uuid}/intermediate" : "${params.output_dir_base}/${params.hisat2_version}/intermediate"
   aligner_validation_dir = (params.ucla_cds_registered_dataset_output) ? "${params.output_dir_base}/${params.hisat2_version}/BAM-${params.hisat2_uuid}/validation" : "${params.output_dir_base}/${params.hisat2_version}/validation"
   aligner_log_dir = "${params.log_output_dir}/process-log/${params.hisat2_version}"
   aligner_qc_dir = (params.ucla_cds_registered_dataset_output) ? "${params.output_dir_base}/${params.hisat2_version}/BAM-${params.hisat2_uuid}/QC" : "${params.output_dir_base}/${params.hisat2_version}/QC"
   take:
      complete_signal //Output bam from previous MarkDuplicatesSpark process to ensure only one Spark process runs at a time
      ich_samples
      ich_samples_validate
      ich_reference_fasta
      ich_reference_index_files
   main:
      input_validation = ich_samples_validate.mix(
            ich_reference_fasta,
            ich_reference_index_files
         )
      validate_input_HISAT2(input_validation)

      validate_input_HISAT2.out.validation_result.collectFile(
         name: 'input_validation.txt',
         storeDir: "${aligner_validation_dir}"
         )
      align_DNA_HISAT2(
         ich_samples,
         ich_reference_fasta,
         ich_reference_index_files.collect()
         )
      run_sort_SAMtools(align_DNA_HISAT2.out.bam, aligner_output_dir, aligner_intermediate_dir, aligner_log_dir)

      remove_intermediate_files(
         run_sort_SAMtools.out.bam_for_deletion,
         "decoy_signal"
         )

      if (!params.mark_duplicates) {
         // It's possible that run_sort_SAMtools may output multiple BAM files which need to be merged
         // only need to merge when !params.mark_duplicates, since  run_MarkDuplicatesSpark_GATK and run_MarkDuplicate_Picard automatically handle multiple BAMs
         run_merge_SAMtools(run_sort_SAMtools.out.bam.collect(), aligner_output_dir, aligner_intermediate_dir, aligner_log_dir)
         och_bam_index = run_merge_SAMtools.out.merged_bam_index
         och_bam = run_merge_SAMtools.out.merged_bam

      } else {
         if (params.enable_spark) {
            //Run MarkduplicatesSpark only after BWA-MEM2 markduplicatesspark completes
            run_MarkDuplicatesSpark_GATK(complete_signal, run_sort_SAMtools.out.bam.collect(), aligner_output_dir, aligner_intermediate_dir, aligner_log_dir, aligner_qc_dir)
            och_bam = run_MarkDuplicatesSpark_GATK.out.bam
            och_bam_index = run_MarkDuplicatesSpark_GATK.out.bam_index
         } else {
            run_MarkDuplicate_Picard(run_sort_SAMtools.out.bam.collect(), aligner_output_dir, aligner_intermediate_dir, aligner_log_dir, aligner_qc_dir)
            och_bam = run_MarkDuplicate_Picard.out.bam
            och_bam_index = run_MarkDuplicate_Picard.out.bam_index
         }
      }
      generate_sha512sum(och_bam_index.mix(och_bam), aligner_output_dir)

      output_validation = och_bam_index

      validate_output_HISAT2(output_validation)

      validate_output_HISAT2.out.validation_result.collectFile(
         name: 'output_validation.txt',
         storeDir: "${aligner_validation_dir}"
         )
   }

