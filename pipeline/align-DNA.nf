#!/usr/bin/env nextflow 

/* 
   TODO:
      - UPDATE DOCKER IMAGE AND COMMAND CALL FOR VALIDATING INPUTS AND OUTPUTS
      - NEEDS LOGGING UPDATE FOR OUTPUS, LOGS AND REPORTS
*/

def docker_image_BWA_and_SAMTools = "blcdsdockerregistry/align-dna:bwa-mem2-2.1_samtools-1.10"
def docker_image_PicardTools = "blcdsdockerregistry/align-dna:picardtools-2.23.3"
def docker_image_sha512sum = "blcdsdockerregistry/align-dna:sha512sum-1.0"
def docker_image_validate_params = "blcdsdockerregistry/align-dna:sha512sum-1.0"

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
      reference_fasta_dict: ${params.reference_fasta_dict}
      reference_fasta_index_files: ${params.reference_fasta_index_files}

   - output: 
      temp_dir: ${params.temp_dir}
      output_dir: ${params.output_dir}
      
   - options:
      save_intermediate_files = ${params.save_intermediate_files}
      cache_intermediate_pipeline_steps = ${params.cache_intermediate_pipeline_steps}
      max_number_of_parallel_jobs = ${params.max_number_of_parallel_jobs}
      bwa_mem_number_of_cpus = ${params.bwa_mem_number_of_cpus}

   Tools Used:
   - BWA-MEM2 and SAMtools: ${docker_image_BWA_and_SAMTools}
   - Picard Tools: ${docker_image_PicardTools}
   - sha512sum: ${docker_image_sha512sum}
   - validate_params: ${docker_image_validate_params}

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
   .into { input_ch_samples_validate; input_ch_samples }

// get the reference and required files for aligning
Channel
   .fromPath(params.reference_fasta)
   .ifEmpty { error "Cannot find reference: ${params.reference_fasta}" }
   .into { input_ch_reference_fasta_validate; input_ch_reference_fasta }

Channel
   .fromPath(params.reference_fasta_dict)
   .ifEmpty { error "Cannot find reference dictionary: ${params.reference_fasta_dict}" }
   .into { input_ch_reference_dict_validate; input_ch_reference_dict }

Channel
   .fromPath(params.reference_fasta_index_files)
   .ifEmpty { error "Cannot find reference index files: ${params.reference_fasta_index_files}" }
   .into { input_ch_reference_index_files_validate; input_ch_reference_index_files }

// get relavent inputs from csv that should be validated
input_ch_samples_validate
   .flatMap { library, lane, read_group_name, read1_fastq, read2_fastq ->
      [read1_fastq, read2_fastq]
   }
   .set { input_ch_2_samples_validate }

process validate_inputs {
   container docker_image_validate_params

   input:
   path(file_to_validate) from input_ch_2_samples_validate.mix(
      input_ch_reference_fasta_validate, 
      input_ch_reference_dict_validate, 
      input_ch_reference_index_files_validate
   )

   output:
      val(true) into output_ch_validate_inputs

   script:
   """
   set -euo pipefail

   #python -m validate -t ${file_to_validate}
   """
}

int number_of_invalid_inputs = output_ch_validate_inputs
         .filter { valid_result ->
            valid_result == false
         }
         .count()
         .get()

// align with bwa mem and convert with samtools
process BWA_mem_SAMTools_Convert_Sam_to_Bam {
   container docker_image_BWA_and_SAMTools
   publishDir path: params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   memory amount_of_memory
   cpus params.bwa_mem_number_of_cpus

   when:
      number_of_invalid_inputs == 0

   // use "each" so the the reference files are passed through for each fastq pair alignment 
   input: 
      tuple(val(library), 
         val(lane),
         val(read_group_name), 
         path(read1_fastq),
         path(read2_fastq) 
      ) from input_ch_samples
      each file(ref_fasta) from input_ch_reference_fasta
      each file(ref_dict) from input_ch_reference_dict
      file(ref_idx_files) from input_ch_reference_index_files.collect()

   // output the lane information in the file name to differentiate bewteen aligments of the same
   // sample but different lanes
   output:
      tuple(val(library), 
         val(lane),
         file("${library}-${lane}.aligned.bam")
      ) into output_ch_BWA_mem_SAMTools_Convert_Sam_to_Bam

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
   container docker_image_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"
   
   publishDir path: params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   memory amount_of_memory
   cpus number_of_cpus

   input:
      tuple(val(library), 
         val(lane),
         file(input_bam)
      ) from output_ch_BWA_mem_SAMTools_Convert_Sam_to_Bam
   
   // the first value of the tuple will be used as a key to group aligned and filtered bams
   // from the same sample and library but different lane together

   // the next steps of the pipeline are merging so using a lane to differentiate between files is no londer needed
   // (files of same lane are merged together) so the lane information is dropped
   output:
      file("${library}-${lane}.sorted.bam") into output_ch_PicardTools_SortSam

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
   container docker_image_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.output_dir, mode: 'copy'

   memory amount_of_memory
   cpus number_of_cpus

   input:
      file(input_bams) from output_ch_PicardTools_SortSam.collect()

   // after marking duplicates, bams will be merged by library so the library name is not needed
   // just the sample name (global variable), do not pass it as a val
   output:
      file("${params.sample_name}.bam") into output_ch_PicardTools_MarkDuplicates

   shell:
   ''' 
   set -euo pipefail

   # add picard option prefix, '--INPUT' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / --INPUT /g' | sed '1s/^/--INPUT /')

   java -Xmx!{params.mem_command_mark_duplicates} -Djava.io.tmpdir=/temp_dir \
      -jar /picard-tools/picard.jar \
      MarkDuplicates \
      --VALIDATION_STRINGENCY LENIENT \
      $INPUT \
      --OUTPUT !{params.sample_name}.bam \
      --METRICS_FILE !{params.sample_name}.mark_dup.metrics \
      --ASSUME_SORT_ORDER coordinate \
      --PROGRAM_RECORD_ID MarkDuplicates
   '''
}

output_ch_PicardTools_MarkDuplicates
   .into { input_ch_PicardTools_BuildBamIndex; input_ch_generate_sha512sum; input_ch_validate_outputs }

// index bams with picard
process PicardTools_BuildBamIndex  {
   container docker_image_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.output_dir, mode: 'copy'

   memory amount_of_memory
   cpus number_of_cpus
   
   input:
      file(input_bam) from input_ch_PicardTools_BuildBamIndex

   // no need for an output channel becuase this is the final stepp
   output:
      file("${input_bam.getName()}.bai") into output_ch_PicardTools_BuildBamIndex

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

output_ch_PicardTools_BuildBamIndex
   .into{ input_ch_2_generate_sha512sum; input_ch_2_validate_outputs }

// produce checksum for the bam and bam index
process generate_sha512sum {    
   container docker_image_sha512sum

   publishDir path: params.output_dir, mode: 'copy'
   
   memory amount_of_memory
   cpus number_of_cpus

   input:
      file(input_file) from input_ch_generate_sha512sum.mix(input_ch_2_generate_sha512sum)

   output:
      file("${input_file.getName()}.sha512sum") into output_ch_generate_sha512sum

   script:
   """
   set -euo pipefail

   sha512sum ${input_file} > ${input_file.getName()}.sha512sum
   """
}

process validate_outputs {
   container docker_image_validate_params

   input:
      path(file_to_validate) from input_ch_validate_outputs
         .mix(output_ch_generate_sha512sum, 
              input_ch_2_validate_outputs,
              Channel.from(params.temp_dir, 
                           params.output_dir))

   script:
   """
   set -euo pipefail

   #python -m validate -t ${file_to_validate}
   """
}
