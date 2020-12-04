#!/usr/bin/env nextflow 

/* 
   TODO:
      - UPDATE DOCKER IMAGE AND COMMAND CALL FOR VALIDATING INPUTS AND OUTPUTS
      - NEEDS LOGGING UPDATE FOR OUTPUS, LOGS AND REPORTS
*/

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

amount_of_memory = amount_of_memory.toString() + " GB"

// Default memory configuration for Picard's Java commands
params.mem_command_sort_sam = "4g"
params.mem_command_mark_duplicates = "4g"
params.mem_command_build_bam_index = "4g"

// cpus for bwa-mem2
// The memory requied by bwa-mem2 increases along with the threads. Here we ensure that each cpu has
// at least 2.5 GB of memory, to avoid out-of-memory failure, unless the number of cpu is defined in
// config.
if (!params.containsKey("bwa_mem_number_of_cpus")) {
   amount_of_memory_int = amount_of_memory.replace(" GB", "") as Integer
   params.bwa_mem_number_of_cpus = (int) Math.min(number_of_cpus, Math.floor(amount_of_memory_int / 2.5))
   if (params.bwa_mem_number_of_cpus < 1) {
      params.bwa_mem_number_of_cpus  = 1
   }
}

// output details of the pipeline run to stdout
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

// get the input fastq pairs
Channel
   .fromPath(params.input_csv)
   .ifEmpty { error "Cannot find input csv: ${params.input_csv}" }
   .splitCsv(header:true)
   .map { row -> 
      def read_group_name = "@RG" +
         "\\tID:" + row.read_group_identifier + ".Seq" + row.lane +
         "\\tCN:" + row.sequencing_center +
         "\\tLB:" + row.library_identifier +
         "\\tPL:" + row.platform_technology +
         "\\tPU:" + row.platform_unit +
         "\\tSM:" + row.sample

      // the library, sample and lane are used as keys downstream to group into 
      // sets of the same key for downstream merging
      return tuple(row.library_identifier,
         row.lane,
         read_group_name,
         row.read1_fastq,
         row.read2_fastq
      )
   }
   .into { ich_samples_validate; ich_samples; ich_samples_output }

// get the reference and required files for aligning
Channel
   .fromPath(params.reference_fasta)
   .ifEmpty { error "Cannot find reference: ${params.reference_fasta}" }
   .into { ich_reference_fasta_validate; ich_reference_fasta }

Channel
   .fromPath(params.reference_fasta_index_files)
   .ifEmpty { error "Cannot find reference index files: ${params.reference_fasta_index_files}" }
   .into { ich_reference_index_files_validate; ich_reference_index_files }

// get relavent inputs from csv that should be validated
ich_samples_validate
   .flatMap { library, lane, read_group_name, read1_fastq, read2_fastq ->
      [read1_fastq, read2_fastq]
   }
   .set { ich_2_samples_validate }

process validate_inputs {
   container dockeri_validate_params

   input:
   path(file_to_validate) from ich_2_samples_validate.mix(
      ich_reference_fasta_validate, 
      ich_reference_index_files_validate
   )

   script:
   """
   set -euo pipefail

   python -m validate -t file-input ${file_to_validate}
   """
}

// align with bwa mem and convert with samtools
process align_BWA_mem_convert_SAM_to_BAM_samtools {
   container dockeri_BWA_and_SAMtools

   publishDir path: params.bam_output_dir,
      enabled: params.save_intermediate_files,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "align_BWA_mem_convert_SAM_to_BAM_samtools/${file(read1_fastq).getName()}/log${file(it).getName()}" }

   memory amount_of_memory
   cpus params.bwa_mem_number_of_cpus

   // use "each" so the the reference files are passed through for each fastq pair alignment 
   input: 
      tuple(val(library), 
         val(lane),
         val(read_group_name), 
         path(read1_fastq),
         path(read2_fastq) 
      ) from ich_samples
      each file(ref_fasta) from ich_reference_fasta
      file(ref_idx_files) from ich_reference_index_files.collect()

   // output the lane information in the file name to differentiate bewteen aligments of the same
   // sample but different lanes
   output:
      tuple(val(library), 
         val(lane),
         file("${library}-${lane}.aligned.bam")
      ) into och_align_BWA_mem_convert_SAM_to_BAM_samtools
      file ".command.*"

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

// sort coordinate order with picard
process PicardTools_SortSam  {
   container dockeri_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"
   
   publishDir path: params.output_dir,
      enabled: params.save_intermediate_files,
      pattern: "*.sorted.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "PicardTools_SortSam/log${file(it).getName()}" }

   memory amount_of_memory
   cpus number_of_cpus

   input:
      tuple(val(library), 
         val(lane),
         file(input_bam)
      ) from och_align_BWA_mem_convert_SAM_to_BAM_samtools
   
   // the first value of the tuple will be used as a key to group aligned and filtered bams
   // from the same sample and library but different lane together

   // the next steps of the pipeline are merging so using a lane to differentiate between files is no londer needed
   // (files of same lane are merged together) so the lane information is dropped
   output:
      file("${library}-${lane}.sorted.bam") into och_PicardTools_SortSam
      file ".command.*"

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
   container dockeri_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.bam_output_dir,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "PicardTools_MarkDuplicates/log${file(it).getName()}" }

   memory amount_of_memory
   cpus number_of_cpus

   input:
      file(input_bams) from och_PicardTools_SortSam.collect()

   // after marking duplicates, bams will be merged by library so the library name is not needed
   // just the sample name (global variable), do not pass it as a val
   output:
      file(bam_output_filename) into och_PicardTools_MarkDuplicates
      file ".command.*"

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

och_PicardTools_MarkDuplicates
   .into { ich_PicardTools_BuildBamIndex; ich_generate_sha512; ich_validate_outputs }

// index bams with picard
process PicardTools_BuildBamIndex  {
   container dockeri_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.bam_output_dir,
      pattern: "*.bam.bai",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "PicardTools_BuildBamIndex/log${file(it).getName()}" }

   memory amount_of_memory
   cpus number_of_cpus
   
   input:
      file(input_bam) from ich_PicardTools_BuildBamIndex

   // no need for an output channel becuase this is the final stepp
   output:
      file("${input_bam.getName()}.bai") into och_PicardTools_BuildBamIndex
      file ".command.*"

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

och_PicardTools_BuildBamIndex
   .into{ ich_2_generate_sha512; ich_2_validate_outputs }

// produce checksum for the bam and bam index
process generate_sha512sum {    
   container dockeri_sha512sum

   publishDir path: params.bam_output_dir, mode: 'copy'
   
   memory amount_of_memory
   cpus number_of_cpus

   input:
      file(input_file) from ich_generate_sha512.mix(ich_2_generate_sha512)

   output:
      file("${input_file.getName()}.sha512") into och_generate_sha512

   script:
   """
   set -euo pipefail

   sha512sum ${input_file} > ${input_file.getName()}.sha512
   """
}

process validate_outputs {
   container dockeri_validate_params

   input:
      path(file_to_validate) from ich_validate_outputs
         .mix(och_generate_sha512, 
              ich_2_validate_outputs,
              Channel.from(params.temp_dir, 
                           params.output_dir))

   script:
   """
   set -euo pipefail

   python -m validate -t file-input ${file_to_validate}
   """
}
