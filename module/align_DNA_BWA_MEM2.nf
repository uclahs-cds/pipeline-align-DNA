// The follwing process runs both alignment and SAM conversion in the same process.
// While this is normally considered to go against the best practices for processes,
// here it actually saves cost, time, and memory to directly pipe the output into 
// samtools due to the large size of the uncompressed SAM files.
include { run_validate; run_validate as validate_output_file } from './validation.nf'
include { run_samtools_sort ; run_index_SAMtools } from './samtools.nf'
include { run_MarkDuplicate_Picard } from './mark_duplicate_picardtools.nf'
include { run_MarkDuplicatesSpark_GATK } from './mark_duplicates_spark.nf'
include { Generate_Sha512sum } from './check_512sum.nf'
include { remove_intermediate_files } from '../external/nextflow-modules/modules/common/intermediate_file_removal/main.nf' addParams(
   options: [
      save_intermediate_files: params.save_intermediate_files,
      output_dir: params.base_output_dir,
      log_output_dir: "${params.log_output_dir}/process-log/${params.bwa_version}"
      ]
   )

process align_DNA_BWA_MEM2 {
   container params.docker_image_bwa_and_samtools
   publishDir path: "${params.base_output_dir}/${params.bwa_version}/intermediate/${task.process.split(':')[1].replace('_', '-')}",
      enabled: params.save_intermediate_files,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: "${params.log_output_dir}/process-log/${params.bwa_version}/${task.process.split(':')[1].replace('_', '-')}",
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
         path("${library}-${lane}.bam"), emit: bam
      path(".command.*")

   script:
   """
   set -euo pipefail

   bwa-mem2 \
      mem \
      -t ${task.cpus} \
      -M \
      -R \"@RG\\tID:${header.read_group_identifier}.Seq${header.lane}\\tCN:${header.sequencing_center}\\tLB:${header.library_identifier}\\tPL:${header.platform_technology}\\tPU:${header.platform_unit}\\tSM:${header.sample}\" \
      ${ref_fasta} \
      ${read1_fastq} \
      ${read2_fastq} | \
   samtools \
      view \
      -@ ${task.cpus} \
      -S \
      -b > \
      ${library}-${lane}.bam
   """
   }

workflow align_DNA_BWA_MEM2_workflow {
   aligner_output_dir = "${params.base_output_dir}/${params.bwa_version}/output"
   aligner_intermediate_dir = "${params.base_output_dir}/${params.bwa_version}/intermediate"
   aligner_validation_dir = "${params.base_output_dir}/${params.bwa_version}/validation"
   aligner_log_dir = "${params.log_output_dir}/process-log/${params.bwa_version}"
   aligner_qc_dir = "${params.base_output_dir}/${params.bwa_version}/QC"
  
   take:
      ich_samples
      ich_samples_validate
      ich_reference_fasta
      ich_reference_index_files
   main:
      run_validate(ich_samples_validate.mix(
         ich_reference_fasta,
         ich_reference_index_files
         ),
         aligner_log_dir
         )

      // change validation file name depending on whether inputs or outputs are being validated
      //val_filename = ${task.process.split(':')[1].replace('_', '-')} == run-validate ? "input_validation.txt" : "output_validation.txt"
      run_validate.out.val_file.collectFile(
         name: 'input_validation.txt',
         storeDir: "${aligner_validation_dir}"
         )
      align_DNA_BWA_MEM2(
         ich_samples,
         ich_reference_fasta,
         ich_reference_index_files.collect()
         )
         
      run_samtools_sort(align_DNA_BWA_MEM2.out.bam, aligner_output_dir, aligner_intermediate_dir, aligner_log_dir)

      remove_intermediate_files(
         run_samtools_sort.out.bam_for_deletion,
         "decoy_signal" 
         ) 
      
      if (!params.mark_duplicates) {
         run_index_SAMtools(run_samtools_sort.out.bam, aligner_output_dir, aligner_intermediate_dir, aligner_log_dir)
         och_bam_index = run_index_SAMtools.out.index
         och_bam = run_samtools_sort.out.bam
      } else {
         if (params.enable_spark) {
            run_MarkDuplicatesSpark_GATK("completion_placeholder", run_samtools_sort.out.bam.collect(), aligner_output_dir, aligner_intermediate_dir, aligner_log_dir, aligner_qc_dir)
            och_bam = run_MarkDuplicatesSpark_GATK.out.bam
            och_bam_index = run_MarkDuplicatesSpark_GATK.out.bam_index
         } else {
            run_MarkDuplicate_Picard(run_samtools_sort.out.bam.collect(), aligner_output_dir, aligner_intermediate_dir, aligner_log_dir, aligner_qc_dir)
            och_bam = run_MarkDuplicate_Picard.out.bam
            och_bam_index = run_MarkDuplicate_Picard.out.bam_index
         }
      }
      Generate_Sha512sum(och_bam_index.mix(och_bam), aligner_output_dir)
      validate_output_file(
         och_bam.mix(
            och_bam_index,
            Channel.from(params.work_dir, params.output_dir)
            ),
            aligner_log_dir
         )
      validate_output_file.out.val_file.collectFile(
         name: 'output_validation.txt',
         storeDir: "${aligner_validation_dir}"
         )
      
      emit:
      complete_signal = och_bam.collect()
   }