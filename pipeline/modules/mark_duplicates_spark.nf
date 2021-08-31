
// mark duplicates with Spark-enabled tool from GATK
process run_MarkDuplicatesSpark_GATK  {
   container params.docker_image_gatk
   containerOptions "--volume ${params.temp_dir}:/temp_dir -u nobody"

   publishDir path: "${bam_output_dir}",
      pattern: "*.bam{,.bai}",
      mode: 'copy'

   publishDir path: "${bam_output_dir}"
      pattern: "*.metrics",
      enabled: params.save_intermediate_files,
      mode: 'copy'

   publishDir path: "${params.log_output_dir}/${task.process.replace(':', '/')}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "log${file(it).getName()}" }

   input:
      path(input_bams)
      val(bam_output_dir)

   // after marking duplicates, bams will be merged by library so the library name is not needed
   // just the sample name (global variable), do not pass it as a val
   output:
      path bam_output_filename, emit: bam
      path bam_index_output_filename, emit: bam_index
      path "${params.sample_name}.mark_dup.metrics"
      path(".command.*")

   beforeScript 'chmod 777 `pwd`'

   shell:
   bam_output_filename = "${params.bam_output_filename}"
   bam_index_output_filename = "${params.bam_output_filename}.bai"
   ''' 
   set -euo pipefail

   # add gatk option prefix, '--input' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / --input /g' | sed '1s/^/--input /')

   gatk MarkDuplicatesSpark \
      --read-validation-stringency LENIENT \
      $INPUT \
      --output !{bam_output_filename} \
      --metrics-file !{params.sample_name}.mark_dup_metrics \
      --program-name MarkDuplicatesSpark \
      --create-output-bam-index \
      --conf 'spark.executor.cores=${task.cpus}' \
      --conf 'spark.local.dir=${spark_temp_dir}' \
      --tmp-dir /temp_dir
   '''
   }
