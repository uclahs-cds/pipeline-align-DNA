
// sort coordinate order with picard
process PicardTools_SortSam  {
   container params.docker_image_picardtools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"
   
   publishDir path: params.output_dir,
      enabled: params.save_intermediate_files,
      pattern: "*.sorted.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "PicardTools_SortSam/log${path(it).getName()}" }

   memory params.amount_of_memory
   cpus params.number_of_cpus

   input:
      tuple(val(library), 
         val(lane),
         path(input_bam)
         )
   
   // the first value of the tuple will be used as a key to group aligned and filtered bams
   // from the same sample and library but different lane together

   // the next steps of the pipeline are merging so using a lane to differentiate between files is no londer needed
   // (files of same lane are merged together) so the lane information is dropped
   output:
      path("${library}-${lane}.sorted.bam")

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

// mark duplicates with picard
process PicardTools_MarkDuplicates  {
   container params.docker_image_picardtools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.bam_output_dir,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "PicardTools_MarkDuplicates/log${path(it).getName()}" }

   memory params.amount_of_memory
   cpus params.number_of_cpus

   input:
      path(input_bams)

   // after marking duplicates, bams will be merged by library so the library name is not needed
   // just the sample name (global variable), do not pass it as a val
   output:
      path(bam_output_filename)

   shell:
   bam_output_filename = params.bam_output_filename
   ''' 
   set -euo pipefail

   # add picard option prefix, '--INPUT' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / --INPUT /g' | sed '1s/^/--INPUT /')

   java -Xmx!{params.mem_command_mark_duplicates} -Djava.io.tmpdir=/temp_dir \
      -jar /picard-tools/picard.jar \
      MarkDuplicates \
      --VALIDATION_STRINGENCY LENIENT \
      $INPUT \
      --OUTPUT !{bam_output_filename} \
      --METRICS_FILE !{params.sample_name}.mark_dup.metrics \
      --ASSUME_SORT_ORDER coordinate \
      --PROGRAM_RECORD_ID MarkDuplicates
   '''
   }

// index bams with picard
process PicardTools_BuildBamIndex  {
   container params.docker_image_picardtools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.bam_output_dir,
      pattern: "*.bam.bai",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "PicardTools_BuildBamIndex/log${path(it).getName()}" }

   memory params.amount_of_memory
   cpus params.number_of_cpus
   
   input:
      path(input_bam)

   // no need for an output channel becuase this is the final stepp
   output:
      path("${input_bam.getName()}.bai")

   script:
   """
   set -euo pipefail

   java -Xmx${params.mem_command_build_bam_index} -Djava.io.tmpdir=/temp_dir \
      -jar /picard-tools/picard.jar \
      BuildBamIndex \
      --VALIDATION_STRINGENCY LENIENT \
      --INPUT ${input_bam} \
      --OUTPUT ${input_bam.getName()}.bai
   """
   }