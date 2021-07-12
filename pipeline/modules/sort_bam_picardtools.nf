
// sort coordinate order with picard
process PicardTools_SortSam  {
   container params.docker_image_picardtools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"
   
   publishDir path: "${intermediate_output_dir}",
      enabled: params.save_intermediate_files,
      pattern: "*.sorted.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "PicardTools_SortSam/${library}-${lane}.log${file(it).getName()}" }

   input:
      tuple(val(library), 
         val(lane),
         path(input_bam)
         )
      val(intermediate_output_dir)
   
   // the first value of the tuple will be used as a key to group aligned and filtered bams
   // from the same sample and library but different lane together

   // the next steps of the pipeline are merging so using a lane to differentiate between files is no longer needed
   // (files of same lane are merged together) so the lane information is dropped
   output:
      path "${library}-${lane}.sorted.bam", emit: bam
      path(".command.*")

   script:
   """
   set -euo pipefail

   java -Xmx${params.mem_command_sort_sam} -Djava.io.tmpdir=/temp_dir \
      -jar /picard-tools/picard.jar \
      SortSam \
      --VALIDATION_STRINGENCY LENIENT \
      --INPUT ${input_bam} \
      --OUTPUT ${library}-${lane}.sorted.bam \
      --SORT_ORDER coordinate
   """
   }
