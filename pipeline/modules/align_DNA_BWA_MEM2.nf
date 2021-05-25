// The follwing process runs both alignment and SAM conversion in the same process.
// While this is normally considered to go against the best practices for processes,
// here it actually saves cost, time, and memory to directly pipe the output into 
// samtools due to the large size of the uncompressed SAM files.

process align_DNA_BWA_MEM2 {
   container params.docker_image_bwa_and_samtools
   publishDir path: params.bam_output_dir,
      enabled: params.save_intermediate_files,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "align_DNA_BWA_MEM2/${file(read1_fastq).getSimpleName()}/${library}-${lane}.log${file(it).getName()}" }

   // use "each" so the the reference files are passed through for each fastq pair alignment 
   input: 
      tuple(val(library),
	 val(header), 
         val(lane),
         val(read_group_name), 
         path(read1_fastq),
         path(read2_fastq) 
         )
      each path(ref_fasta)
      path(ich_reference_index_files)

   // output the lane information in the file name to differentiate bewteen aligments of the same
   // sample but different lanes
   output:
      tuple val(library), 
         val(lane),
         path("${library}-${lane}.bam"), emit: bam
      path(".command.*")

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
      ${library}-${lane}.bam
   """
   }
