
// mark duplicates with picard
process run_MarkDuplicatesSpark_GATK  {
   container 'mygatk:1.0.0'
   containerOptions "--volume ${params.temp_dir}:/temp_dir -u sparkuser"

   publishDir path: "${bam_output_dir}",
      pattern: "*.bam{,.bai}",
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
      --tmp-dir /temp_dir

   '''
   }
