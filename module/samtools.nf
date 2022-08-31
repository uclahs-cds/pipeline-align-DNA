include { generate_standard_filename } from '../external/nextflow-modules/modules/common/generate_standardized_filename/main.nf'

// sort coordinate or queryname order with samtools
process run_sort_SAMtools  {
   container params.docker_image_samtools
   
   publishDir path: "${intermediate_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      enabled: params.save_intermediate_files && params.mark_duplicates,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: "${log_output_dir}/${task.process.split(':')[-1].replace('_', '-')}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${library}/${lane}/log${file(it).getName()}" }

   input:
      tuple(val(library), 
         val(lane),
         path(input_bam)
         )
      val(bam_output_dir)
      val(intermediate_output_dir)
      val(log_output_dir)
   
   /** the first value of the tuple will be used as a key to group aligned and filtered bams
   * from the same sample and library but different lane together

   * the next steps of the pipeline are merging so using a lane to differentiate between files is no longer needed
   * (files of same lane are merged together) so the lane information is dropped
   */
   
   output:
      path "${bam_output_filename}", emit: bam
      path input_bam, emit: bam_for_deletion
      path(".command.*")

   script:

   bam_output_filename = params.bam_output_filename.replaceAll('.bam$', "${library}-${lane}.sorted.bam")
   
   /** 
   Determine sort order based on markduplicates process: queryname for spark and coordinate for Picard

   If params.mark_duplicates=false, then sort_order="", thus samtools sort will sort by coordinates
   If params.mark_duplicates=true, then sort by queryname if params.enable_spark=true, since this is what MarkDuplicatesSpark expects.
      otherwise sort by coordinates since this is what Picard MarkDuplicates expects.
   */

   sort_order = (params.mark_duplicates && params.enable_spark) ? "-n" : ""

   if ("-n" == sort_order) {
      println("Sorting by queryname for MarkDuplicatesSpark (sort_order=${sort_order})")
   } else {
      println("Sorting by coordinates (sort_order=${sort_order})")
   }
   
   """
   set -euo pipefail

   samtools sort \
    -@ ${task.cpus} \
    -O bam \
    -o ${bam_output_filename} \
    ${sort_order} \
    ${input_bam}
   """
   }

/* When params.mark_duplicates=false, it's possible that run_sort_SAMtools may output multiple BAM files which 
need to be merged.  When params.mark_duplicates=true, merging is not needed because run_MarkDuplicatesSpark_GATK 
and run_MarkDuplicate_Picard automatically handle multiple BAMs.
*/
process run_merge_SAMtools  {
   container params.docker_image_samtools
   
   publishDir path: "${bam_output_dir}",
      pattern: "${merged_bam_output_filename}{,.bai}",
      mode: 'copy'

   publishDir path: "${log_output_dir}/${task.process.split(':')[-1].replace('_', '-')}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "log${file(it).getName()}" }

   input:
      // outputs from run_sort_SAMtools
      path(bam) // bam file(s)
      val(bam_output_dir) //directory of bam
      val(intermediate_output_dir)
      val(log_output_dir)
   
   output:
      path "${merged_bam_output_filename}", emit: merged_bam
      path "${merged_bam_output_filename}.bai", emit: merged_bam_index
      path(".command.*")

   script:

   merged_bam_output_filename = "${params.bam_output_filename}"

   """
   set -euo pipefail

   samtools merge \
    --threads ${task.cpus} \
    --write-index \
    -o ${merged_bam_output_filename}##idx##${merged_bam_output_filename}.bai \
    ${bam}
   """
   }

