
// sort coordinate or queryname order with picard
process run_SortSam_Picard  {
   container params.docker_image_picardtools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"
   
   publishDir path: "${intermediate_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      enabled: params.save_intermediate_files && params.mark_duplicates,
      pattern: "*.{bam,bai}",
      mode: 'copy',
      saveAs: { filename -> (file(filename).getExtension() == "bai") ? "${file(filename).baseName}.bam.bai" : "${filename}" }

   publishDir path: "${bam_output_dir}",
      enabled: !params.mark_duplicates,
      pattern: "*.{bam,bai}",
      mode: 'copy',
      saveAs: { filename -> (file(filename).getExtension() == "bai") ? "${file(filename).baseName}.bam.bai" : "${filename}" }

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
      path "*.bam", emit: bam
      path "*.bai", emit: bam_index optional true
      path(".command.*")

   script:
   /** 
   * Determine sort order based on markduplicates process: queryname for spark and coordinate for Picard
   * Determine filename based on whether markduplicates processes are enabled.
   * Index output file if sorting is the final step in the pipeline (if markduplicates disabled)
   */

   if (!params.mark_duplicates) {
         sort_order = "coordinate"
         bam_output_filename = "${params.bam_output_filename}"
         index = true
      } else {
         sort_order = (params.enable_spark) ? "queryname" : "coordinate"
         bam_output_filename = "${library}-${lane}.sorted.bam"
         index = false
      }
   
   """
   set -euo pipefail

   java -Xmx${task.memory.getMega()}m -Djava.io.tmpdir=/temp_dir \
      -jar /usr/local/share/picard-slim-2.26.10-0/picard.jar \
      SortSam \
      --VALIDATION_STRINGENCY LENIENT \
      --INPUT ${input_bam} \
      --OUTPUT ${bam_output_filename} \
      --SORT_ORDER ${sort_order} \
      --CREATE_INDEX ${index}
   """
   }
