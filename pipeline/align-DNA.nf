#!/usr/bin/env nextflow 

/* 
   TODO:
      - NEEDS LOGGING UPDATE FOR OUTPUS, LOGS AND REPORTS
      - NEEDS A VALIDATION PROCESS/STEP
      - AFTER SEP 2015 PICARD TOOLS MARK DUPS IS LIBRARY AWARE
      AS A RESULT, WE CAN COMBINE MERGE AND MARK DUPS STEPS INTO 
      1 PROCESS
*/

def docker_image_BWA_and_SAMTools = "blcdsdockerregistry/align-dna:bwa-0.7.17_samtools-1.10"
def docker_image_PicardTools = "blcdsdockerregistry/align-dna:picardtools-2.23.3"
def docker_image_sha512sum = "blcdsdockerregistry/align-dna:sha512sum-1.0"
def docker_image_validate_params = "blcdsdockerregistry/align-dna:sha512sum-1.0"

// resource information
def number_of_cpus = Runtime.getRuntime().availableProcessors()
def total_memory = ((int) ((java.lang.management.ManagementFactory.getOperatingSystemMXBean()
   .getTotalPhysicalMemorySize() / (1024.0 * 1024.0 * 1024.0)) * 0.9)).toString() + " GB"

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

   Tools Used:
   - BWA and SAMtools: ${docker_image_BWA_and_SAMTools}
   - Picard Tools: ${docker_image_PicardTools}
   - sha512sum: ${docker_image_sha512sum}
   - validate_params: ${docker_image_validate_params}
   ------------------------------------
   Executing workflow...
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

int number_of_invalid_inputs = output_ch_validate_inputs.collect()
         .filter { valid_result ->
            valid_result == false
         }
         .count()
         .get()

// align with bwa mem and convert with samtools
process BWA_mem_SAMTools_Convert_Sam_to_Bam {
   container docker_image_BWA_and_SAMTools

   publishDir path: params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   when:
      number_of_invalid_inputs == 0

   memory total_memory
   cpus number_of_cpus

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

   bwa \
      mem \
      -t ${number_of_cpus} \
      -M \
      -R "${read_group_name}" \
      ${ref_fasta} \
      ${read1_fastq} \
      ${read2_fastq} | \
   samtools \
      view \
      -@ ${number_of_cpus} \
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

   memory total_memory
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
      tuple(val(library), 
         file("${library}-${lane}.sorted.bam")
      ) into output_ch_PicardTools_SortSam

   script:
   """
   set -euo pipefail

   java -Xmx4g -Djava.io.tmpdir=/temp_dir \
      -jar /picard-tools/picard.jar \
      SortSam \
      --VALIDATION_STRINGENCY LENIENT \
      --INPUT ${input_bam} \
      --OUTPUT ${library}-${lane}.sorted.bam \
      --SORT_ORDER coordinate
   """
}

// group the aligned and filtered bams by the same library
// and send those outputs to be merged if there is more than one bam per group (this means there are multiple lanes)
output_ch_PicardTools_SortSam
	.groupTuple()
   .branch { library, output_SortSam_bams ->
		input_ch_PicardTools_MergeSamFiles_across_lanes: output_SortSam_bams.size() > 1

      // when grouping in values besides the first value passed, the key become wrapped in an additional tuple
      // and b/c we know that the library and bams are a tuple  of size 1 and the downstream input requires a file and library
      // that are not wrapped in a tuple, we just get the 1st element of the tuple, the files themselves

      // samples that are only aligned to one lane 
		input_ch_PicardTools_MergeSamFiles_across_libraries_PicardTools_MarkDuplicates: output_SortSam_bams.size() <= 1
         return tuple(library, output_SortSam_bams.get(0))
	}
	.set { output_ch_2_PicardTools_SortSam }

// merge bams from across lanes from the same library with picard
process PicardTools_MergeSamFiles_across_lanes  {
   container docker_image_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"
   
   publishDir path: params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   memory total_memory
   cpus number_of_cpus

   input:
      tuple(val(library),
         file(input_bams)
      ) from output_ch_2_PicardTools_SortSam.input_ch_PicardTools_MergeSamFiles_across_lanes

   output:
      tuple(val(library),
         file("${library}.merged-lanes.bam")
      ) into output_ch_PicardTools_MergeSamFiles_across_lanes

   shell:
   '''
   set -euo pipefail

   # add picard option prefix, 'INPUT=' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / --INPUT /g' | sed '1s/^/--INPUT /')

   java -Xmx4g -Djava.io.tmpdir=/temp_dir \
      -jar /picard-tools/picard.jar \
      MergeSamFiles \
      --USE_THREADING true \
      --VALIDATION_STRINGENCY LENIENT \
      $INPUT \
      --OUTPUT !{library}.merged-lanes.bam
   '''
}

output_ch_PicardTools_MergeSamFiles_across_lanes
   .mix(output_ch_2_PicardTools_SortSam.input_ch_PicardTools_MergeSamFiles_across_libraries_PicardTools_MarkDuplicates)
   .set { input_ch_PicardTools_MarkDuplicates }

// mark duplicates with picard
process PicardTools_MarkDuplicates  {
   container docker_image_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   memory total_memory
   cpus number_of_cpus

   input:
      tuple(val(library), 
         file(input_bam)
      ) from input_ch_PicardTools_MarkDuplicates

   // after marking duplicates, bams will be merged by library so the library name is not needed
   // just the sample name (global variable), do not pass it as a val
   output:
      file("${library}.mark_dup.bam") into output_ch_PicardTools_MarkDuplicates
      file("${library}.mark_dup.metrics")

   shell:
   """
   set -euo pipefail

   java -Xmx4g -Djava.io.tmpdir=/temp_dir/ \
      -jar /picard-tools/picard.jar \
      MarkDuplicates \
      --VALIDATION_STRINGENCY LENIENT \
      --INPUT ${input_bam} \
      --OUTPUT ${library}.mark_dup.bam \
      --METRICS_FILE ${library}.mark_dup.metrics \
      --PROGRAM_RECORD_ID MarkDuplicates
   """
}

// copy into 2 channels and use one to get the size of the channel
output_ch_PicardTools_MarkDuplicates
   .into { output_ch_PicardTools_MarkDuplicates_count; output_ch_2_PicardTools_MarkDuplicates }

//  the number of mark duplicate bams == the number of libraries 
int num_of_libraries = output_ch_PicardTools_MarkDuplicates_count
	.count()
	.get()

// send to merge across lanes if more than one library
output_ch_2_PicardTools_MarkDuplicates
   .branch { bams ->
		input_ch_PicardTools_MergeSamFiles_across_libraries: num_of_libraries > 1
      
      // if only one library
		input_ch_PicardTools_BuildBamIndex: true
	}
	.set { output_ch_3_PicardTools_MarkDuplicates }

// merge bams from the same library with picard
process PicardTools_MergeSamFiles_across_libraries  {
   container docker_image_PicardTools
   containerOptions "--volume ${params.temp_dir}:/temp_dir"

   publishDir path: params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   memory total_memory
   cpus number_of_cpus

   input:
      file(input_bams) from output_ch_3_PicardTools_MarkDuplicates.input_ch_PicardTools_MergeSamFiles_across_libraries.collect()

   output:
      file("${params.sample_name}.merged-libraries.bam") into output_ch_PicardTools_MergeSamFiles_across_libraries

   shell:
   '''
   set -euo pipefail

   # add picard option prefix, 'INPUT=' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / --INPUT /g' | sed '1s/^/--INPUT /')

   java -Xmx4g -Djava.io.tmpdir=/temp_dir \
      -jar /picard-tools/picard.jar \
      MergeSamFiles \
      --USE_THREADING true \
      --VALIDATION_STRINGENCY LENIENT \
      $INPUT \
      --OUTPUT !{params.sample_name}.merged-libraries.bam
   '''
}

output_ch_3_PicardTools_MarkDuplicates.input_ch_PicardTools_BuildBamIndex
   .mix(output_ch_PicardTools_MergeSamFiles_across_libraries)
   .into { input_ch_2_PicardTools_BuildBamIndex; input_ch_generate_sha512sum;  }

// index bams with picard
process PicardTools_BuildBamIndex  {
   container docker_image_PicardTools
   containerOptions "--volume ${params.java_temp_dir}:/java_temp_dir"

   publishDir path: params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   memory total_memory
   cpus number_of_cpus
   
   input:
      file(input_bam) from input_ch_2_PicardTools_BuildBamIndex

   // no need for an output channel becuase this is the final stepp
   output:
      file("${input_bam.getName()}.bai") into output_ch_PicardTools_BuildBamIndex

   script:
   """
   set -euo pipefail

   java -Xmx4g -Djava.io.tmpdir=/java_temp_dir \
      -jar /picard-tools/picard.jar \
      BuildBamIndex \
      --VALIDATION_STRINGENCY LENIENT \
      --INPUT ${input_bam} \
      --OUTPUT ${input_bam.getName()}.bai
   """
}

output_ch_PicardTools_BuildBamIndex
   .mix(input_ch_generate_sha512sum)
   .into{ input_ch_2_generate_sha512sum; input_ch_validate_outputs }

// produce checksum for the bam and bam index
process generate_sha512sum {    
   container docker_image_sha512sum

   publishDir path: params.output_dir, mode: 'copy'
   
   memory total_memory
   cpus number_of_cpus

   input:
      file(input_file) from input_ch_2_generate_sha512sum

   output:
      tuple(file("${input_file}"),
        file("${input_file.getName()}.sha512sum")) into output_ch_generate_sha512sum

   shell:
   """
   set -euo pipefail

   sha512sum ${input_file} > ${input_file.getName()}.sha512sum
   """
}

// output validation
output_ch_generate_sha512sum
   .flatten()
   .mix(input_ch_validate_outputs)
   .set { input_ch_2_validate_outputs }

process validate_outputs {
   container docker_image_validate_params

   input:
      path(file_to_validate) from input_ch_2_validate_outputs

   script:
   """
   set -euo pipefail

   #python -m validate -t ${file_to_validate}
   """
}