
process Align_BWA_mem_convert_SAM_to_BAM_samtools {
   container params.docker_image_bwa_and_samtools
   publishDir path: params.bam_output_dir,
      enabled: params.save_intermediate_files,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "align_BWA_mem_convert_SAM_to_BAM_samtools/${file(read1_fastq).getSimpleName()}/log${file(it).getName()}" }

   memory params.amount_of_memory
   cpus params.bwa_mem_number_of_cpus

   // use "each" so the the reference files are passed through for each fastq pair alignment 
   input: 
      tuple(val(library), 
         val(lane),
         val(read_group_name), 
         path(read1_fastq),
         path(read2_fastq) 
      )
      each file(ref_fasta)
      file(ich_reference_index_files)

   // output the lane information in the file name to differentiate bewteen aligments of the same
   // sample but different lanes
   output:
      tuple(val(library), 
         val(lane),
         file("${library}-${lane}.aligned.bam")
      )
      //file ".command.*"
      // into och_align_BWA_mem_convert_SAM_to_BAM_samtools

   script:
   """
   set -euo pipefail

   bwa-mem2 \
      mem \
      -t ${task.cpus} \
      -M \
      -R "${read_group_name}" \
      ${ref_fasta} \
      ${read1_fastq} \
      ${read2_fastq} | \
   samtools \
      view \
      -@ ${task.cpus} \
      -S \
      -b > \
      ${library}-${lane}.aligned.bam
   """
}