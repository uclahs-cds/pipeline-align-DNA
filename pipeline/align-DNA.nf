
nextflow.enable.dsl=2

// TODO: Remove this later?
def dockeri_BWA_and_SAMtools = "blcdsdockerregistry/align-dna:${params.bwa_version}_samtools-1.10"
def dockeri_PicardTools = "blcdsdockerregistry/align-dna:picardtools-2.23.3"
def dockeri_sha512sum = "blcdsdockerregistry/align-dna:sha512sum-1.0"
def dockeri_validate_params = "blcdsdockerregistry/validate:1.0.0"

// resource information
def number_of_cpus = (int) (Runtime.getRuntime().availableProcessors() / params.max_number_of_parallel_jobs)
if (number_of_cpus < 1) {
   number_of_cpus = 1
}
def amount_of_memory = ((int) (((java.lang.management.ManagementFactory.getOperatingSystemMXBean()
   .getTotalPhysicalMemorySize() / (1024.0 * 1024.0 * 1024.0)) * 0.9) / params.max_number_of_parallel_jobs ))
if (amount_of_memory < 1) {
   amount_of_memory = 1
}

params.amount_of_memory = amount_of_memory.toString() + " GB"

// Default memory configuration for Picard's Java commands
params.mem_command_sort_sam = "4g"
params.mem_command_mark_duplicates = "4g"
params.mem_command_build_bam_index = "4g"

// cpus for bwa-mem2
// The memory requied by bwa-mem2 increases along with the threads. Here we ensure that each cpu has
// at least 2.5 GB of memory, to avoid out-of-memory failure, unless the number of cpu is defined in
// config.
if (!params.containsKey("bwa_mem_number_of_cpus")) {
   params.amount_of_memory_int = amount_of_memory.replace(" GB", "") as Integer
   params.bwa_mem_number_of_cpus = (int) Math.min(number_of_cpus, Math.floor(amount_of_memory_int / 2.5))
   if (params.bwa_mem_number_of_cpus < 1) {
      params.bwa_mem_number_of_cpus  = 1
   }
}


log.info """\
   ===================================
   P I P E L I N E - A L I G N - D N A
   ===================================
   Boutros Lab

   Current Configuration:
   - input: 
      sample_name: ${params.sample_name}
      input_csv: ${params.input_csv}
      reference_fasta: ${params.reference_fasta}
      reference_fasta_index_files: ${params.reference_fasta_index_files}

   - output: 
      temp_dir: ${params.temp_dir}
      output_dir: ${params.output_dir}
      bam_output_dir: ${params.bam_output_dir}
      bam_output_filename: ${params.bam_output_filename}
      log_output_dir: ${params.log_output_dir}
      
   - options:
      save_intermediate_files = ${params.save_intermediate_files}
      cache_intermediate_pipeline_steps = ${params.cache_intermediate_pipeline_steps}
      max_number_of_parallel_jobs = ${params.max_number_of_parallel_jobs}
      bwa_mem_number_of_cpus = ${params.bwa_mem_number_of_cpus}
      blcds_registered_dataset_input = ${params.blcds_registered_dataset_input}
      blcds_registered_dataset_output = ${params.blcds_registered_dataset_output}

   Tools Used:
   - BWA-MEM2 and SAMtools: ${dockeri_BWA_and_SAMtools}
   - Picard Tools: ${dockeri_PicardTools}
   - sha512sum: ${dockeri_sha512sum}
   - validate_params: ${dockeri_validate_params}

   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

include { validate_ichsamples } from './modules/validation-ichsamples.nf'
include { aligndna } from './modules/aligndna-processes.nf'// addParams(bam_output_dir: params.bam_output_dir, log_output_dir: params.log_output_dir)

log.info """\
p1: ${params.bam_output_dir}
p2: ${params.log_output_dir}
"""

workflow {
    ch_bwamem_files = validate_ichsamples()
    //ch_bwamem_files = channel.fromPath([ich_samples, ich_reference_fasta, ich_reference_index_files.collect()])
    aligndna(ch_bwamem_files)
}