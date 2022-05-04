/*
   Nextflow module for marking duplicates using Spark. 
   The containerOptions specifying user 'nobody' allow for Spark to be run without root access.
   The beforeScript command allows the user 'nobody' to create the output files into the working directory.

   Input:
   completion_signal: output bam from previous markduplicatesspark process to ensure 
      only one spark process runs at a time
*/
process run_MarkDuplicatesSpark_GATK  {
   container params.docker_image_gatk
   containerOptions "--volume ${params.work_dir}:/temp_dir --volume ${params.spark_temp_dir}:/spark_temp_dir -u nobody"

   publishDir path: "${bam_output_dir}",
      pattern: "*.bam{,.bai}",
      mode: 'copy'

   publishDir path: "${qc_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      pattern: "*.metrics",
      enabled: params.save_intermediate_files,
      mode: 'copy'

   publishDir path: "${log_output_dir}/${task.process.split(':')[1].replace('_', '-')}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "log${file(it).getName()}" }

   input:
      val(completion_signal)
      path(input_bams)
      val(bam_output_dir)
      val(intermediate_output_dir)
      val(log_output_dir)
      val(qc_output_dir)

   // after marking duplicates, bams will be merged by library so the library name is not needed
   // just the sample name (global variable), do not pass it as a val
   output:
      path bam_output_filename, emit: bam
      path "*.bai", emit: bam_index
      path "${params.sample_name}.mark_dup.metrics"
      path(".command.*")

   //Update tempdir permissions for user 'nobody'
   beforeScript "chmod 777 `pwd`; \
      if [[ ! -d ${params.work_dir} ]]; \
      then \
         mkdir -p ${params.work_dir}; \
         chmod 777 ${params.work_dir}; \
      else \
         if [[ ! `stat -c %a ${params.work_dir}` == 777 ]]; \
         then \
            chmod 777 ${params.work_dir}; \
         fi; \
      fi; \
      if [[ ! -d ${params.spark_temp_dir} ]]; \
      then \
         mkdir -p ${params.spark_temp_dir}; \
         chmod 777 ${params.spark_temp_dir}; \
      else \
         if [[ ! `stat -c %a ${params.spark_temp_dir}` == 777 ]]; \
         then \
            chmod 777 ${params.spark_temp_dir}; \
         fi; \
      fi"

   shell:
   bam_output_filename = "${params.bam_output_filename}"
   ''' 
   set -euo pipefail

   # add gatk option prefix, '--input' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / --input /g' | sed '1s/^/--input /')

   gatk --java-options "-Djava.io.tmpdir=/temp_dir" \
      MarkDuplicatesSpark \
      --read-validation-stringency LENIENT \
      $INPUT \
      --output !{bam_output_filename} \
      --metrics-file !{params.sample_name}.mark_dup.metrics \
      --program-name MarkDuplicatesSpark \
      --create-output-bam-index \
      --conf 'spark.executor.cores=${task.cpus}' \
      --conf 'spark.local.dir=/spark_temp_dir' \
      --tmp-dir /temp_dir
   '''
   }
