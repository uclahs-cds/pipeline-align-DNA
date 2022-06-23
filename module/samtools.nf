// sort coordinate or queryname order with samtools
process run_sort_SAMtools  {
   container params.docker_image_samtools
   
   publishDir path: "${intermediate_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      enabled: params.save_intermediate_files && params.mark_duplicates,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: "${bam_output_dir}",
      enabled: !params.mark_duplicates,
      pattern: "${bam_output_filename}",
      mode: 'copy',
      saveAs: { filename -> "${filename}" }

   publishDir path: "${log_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
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
   /** 
   * Determine sort order based on markduplicates process: queryname for spark and coordinate for Picard
   * Determine filename based on whether markduplicates processes are enabled.
   * Index output file if sorting is the final step in the pipeline (if markduplicates disabled)
   */

   if (!params.mark_duplicates) {
         sort_order = "" /** Empty for sorting by coordinate*/
         bam_output_filename = "${library}-${lane}.sorted.bam"
      } else { 
         sort_order = (params.enable_spark) ? "-n" : "" /** -n to sort by qname*/
         bam_output_filename = "${library}-${lane}.sorted.bam"
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

process run_index_SAMtools  {
   container params.docker_image_samtools
   
   publishDir path: "${intermediate_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      enabled: params.save_intermediate_files,
      pattern: "*.bai",
      mode: 'copy' 

    publishDir path: "${bam_output_dir}",
      enabled: !params.mark_duplicates,
      pattern: "*.bai",
      mode: 'copy' 

   publishDir path: "${log_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "log${file(it).getName()}" }

   input:
      path(bam)
      val(bam_output_dir)
      val(intermediate_output_dir)
      val(log_output_dir)
   
   output:
      path "*.bai", emit: index
      path(".command.*")

   script:
   """
   set -euo pipefail

   samtools index \
    -@ ${task.cpus} \
    ${bam}
   """
   }


   process run_merge_SAMtools  {
   container params.docker_image_samtools
   
   publishDir path: "${intermediate_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      // do we need `&& params.mark_duplicates` ?
      enabled: params.save_intermediate_files && params.mark_duplicates,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: "${log_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${library}/${lane}/log${file(it).getName()}" }

   input:
      // outputs from run_sort_SAMtools
      path(bam) // bam file(s)
      val(bam_output_dir) //directory of bam
      val(intermediate_output_dir)
      val(log_output_dir)
   
   output:
      path "*.bam", emit: merged_bam
      path(".command.*")

   script:

   merged_bam_output_filename = "${library}-sorted-merged.bam"

   """
   set -euo pipefail

   samtools merge \
    --threads ${task.cpus} \
    -o ${merged_bam_output_filename} \
    ${bam}
   """
   }

