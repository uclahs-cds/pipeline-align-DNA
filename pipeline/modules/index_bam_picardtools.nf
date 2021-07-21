
// index bams with picard
process run_BuildBamIndex_Picard  {
   container params.docker_image_picardtools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: "${bam_output_dir}",
      pattern: "*.bam.bai",
      mode: 'copy'

   publishDir path: "${bam_log_output_dir}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "run_BuildBamIndex_Picard/log${file(it).getName()}" }
   
   input:
      path(input_bam)
      val(bam_output_dir)
      val(bam_log_output_dir)

   // no need for an output channel becuase this is the final stepp
   output:
      path "${input_bam.getName()}.bai", emit: bai
      path(".command.*")

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
