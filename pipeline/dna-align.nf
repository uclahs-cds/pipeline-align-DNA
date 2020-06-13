#!/usr/bin/env nextflow 

/* 
   TODO:
      - NEEDS TESTING FOR ALIGNING MULTIPLE BAMS IN PARALLEL
      - NEEDS SLURM OPTIMIZATION AND CONFIGURATION
*/
def docker_image_BWA = "blcdsdockerregistry/bwa:0.7.15"
def docker_image_SAMTools = "blcdsdockerregistry/samtools:1.3"
def docker_image_PicardTools = "blcdsdockerregistry/picard-tools:1.130"

// output details of the pipeline run to stdout
log.info """\
   =============================================
   C A L L - D N A - A L I G N  P I P E L I N E
   =============================================
   Boutros Lab

   Current Configuration:
   - input: 
      input_csv: ${params.input_csv}
      reference_fasta: ${params.reference_fasta}
      reference_fasta_dict: ${params.reference_fasta_dict}
      reference_fasta_index_files: ${params.reference_fasta_index_files}

   - output: 
      output_dir: ${params.output_dir}
      
   - options:
      save_intermediate_files = ${params.save_intermediate_files}

   Tools Used:
   - BWA: ${docker_image_BWA}
   - Samtools: ${docker_image_SAMTools}
   - Picard Tools: ${docker_image_PicardTools}
   
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
         "\tID:" + row.read_group_identifier + ".Seq" + row.lane +
         "\tCN:" + row.sequencing_center +
         "\tLB:" + row.library_identifier +
         "\tPL:" + row.platfrom_technology +
         "\tPU:" + row.platform_unit +
         "\tSM:" + row.sample

      // the library, sample and lane are used as keys downstream to group into 
      // sets of the same key for downstream merging
      return tuple(row.library_identifier,
         row.sample , 
         row.lane,
         read_group_name,
         row.read1_fastq,
         row.read2_fastq
      )
   }
   .set { input_ch_samples }

// get the reference and required files for aligning
Channel
   .fromPath(params.reference_fasta)
   .ifEmpty { error "Cannot find reference: ${params.reference_fasta}" }
   .set { input_ch_reference_fasta }

Channel
   .fromPath(params.reference_fasta_dict)
   .ifEmpty { error "Cannot find reference dictionary: ${params.reference_fasta_dict}" }
   .set { input_ch_reference_dict }

Channel
   .fromPath(params.reference_fasta_index_files)
   .ifEmpty { error "Cannot find reference index files: ${params.reference_fasta_index_files}" }
   .set { input_ch_reference_index_files }
   
// align with bwa mem
process BWA_mem {
   container docker_image_BWA

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   // use "each" so the the reference files are passed through for each fastq pair alignment 
   input: 
      tuple(val(library), 
         val(sample),
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
         val(sample),
         val(lane),
         file("${library}-${sample}-${lane}.aligned.sam")
      ) into output_ch_bwa_mem

   script:
   """
   set -euo pipefail

   bwa \
      mem \
      -t 32 \
      -M \
      -R "${read_group_name}" \
      ${ref_fasta} \
      ${read1_fastq} \
      ${read2_fastq} > \
      ${library}-${sample}-${lane}.aligned.sam
   """
}

// convert with samtools
process SAMTools_convert_sam_to_bam {
   container docker_image_SAMTools

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input: 
      tuple(val(library), 
         val(sample), 
         val(lane),
         file(input_sam)
      ) from output_ch_bwa_mem

   output:
      tuple (val(library), 
         val(sample), 
         val(lane), 
         file("${library}-${sample}-${lane}.aligned.bam")
      ) into output_ch_SAMTools_convert_sam_to_bam

   script:
   """
   set -euo pipefail

   samtools \
      view \
      -@ 4 \
      -S \
      -b \
      ${input_sam} > \
      ${library}-${sample}-${lane}.aligned.bam
   """
}

// sort coordinate order with picard
process PicardTools_SortSam  {
   container docker_image_PicardTools
   
   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(library), 
         val(sample),
         val(lane),
         file(input_bam)
      ) from output_ch_SAMTools_convert_sam_to_bam
   
   output:
      tuple(val(library), 
         val(sample),
         val(lane),
         file("${library}-${sample}-${lane}.sorted.bam")
      ) into output_ch_PicardTools_SortSam

   script:
   """
   set -euo pipefail

   java -Xmx14g -jar /picard-tools/picard.jar \
      SortSam \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${library}-${sample}-${lane}.sorted.bam \
      SORT_ORDER=coordinate
   """
}

// mark duplicates with picard
process PicardTools_MarkDuplicates  {
   container docker_image_PicardTools

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(library), 
         val(sample), 
         val(lane),
         file(input_bam)
      ) from output_ch_PicardTools_SortSam

   // the first value of the tuple will be used as a key to group aligned and filtered bams
   // from the same sample and library but different lane together

   // the next steps of the pipeline are merging so using a lane to differentiate between files is no londer needed
   // (files of same lane are merged together) so the lane information is dropped
   output:
      tuple(val("${library}-${sample}"),
         val(library), 
         file("${library}-${sample}-${lane}.mark_dup.bam")
      ) into output_ch_PicardTools_MarkDuplicates
      file("${library}-${sample}-${lane}.mark_dup.metrics")

   script:
   """
   set -euo pipefail

   java -Xmx10g -jar /picard-tools/picard.jar \
      MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${library}-${sample}-${lane}.mark_dup.bam \
      METRICS_FILE=${library}-${sample}-${lane}.mark_dup.metrics \
      PROGRAM_RECORD_ID=MarkDuplicates
   """
}

// group the aligned and filtered bams the same sample and library
// and send those outputs to be merged if there is more than one bam per group
output_ch_PicardTools_MarkDuplicates
	.groupTuple()
   .branch { sample_and_library_name, library, mark_dups_bams ->
		input_ch_PicardTools_MergeSamFiles_from_same_sample_and_library: mark_dups_bams.size() > 1

      // when grouping in values besides the first value passed, the key become wrapped in an additional tuple
      // and b/c we know that the library and bams are a tuple  of size 1 and the downstream input requires a file and library
      // that are not wrapped in a tuple, we just get the 1st element of the tuple, the file and library themselves

      // samples that are only aligned to one lane will be merged with other samples from the same library
      //  or indexed if there are no other samples from the same library
		input_ch_PicardTools_MergeSamFiles_from_same_library_PicardTools_BuildBamIndex: mark_dups_bams.size() <= 1
         return tuple(library.get(0), mark_dups_bams.get(0))
	}
	.set { output_ch_2_PicardTools_MarkDuplicates }

// merge bams from the same library and sample with picard
process PicardTools_MergeSamFiles_from_same_sample_and_library  {
   container docker_image_PicardTools
   
   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(sample_and_library), 
         val(library),
         file(input_bams)
      ) from output_ch_2_PicardTools_MarkDuplicates.input_ch_PicardTools_MergeSamFiles_from_same_sample_and_library

   // the next steps of the pipeline are merging so using a sample to differentiate between files is no londer needed
   // (files of same sample are merged together) so the sample information is dropped
   output:
      tuple(val(library),
         file("${sample_and_library}.merged.bam")
      ) into output_ch_PicardTools_MergeSamFiles_from_same_sample_and_library

   shell:
   '''
   set -euo pipefail

   # add picard option prefix, 'INPUT=' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / INPUT=/g' | sed '1s/^/INPUT=/')

   java -Xmx10g -jar /picard-tools/picard.jar \
      MergeSamFiles \
      USE_THREADING=true \
      VALIDATION_STRINGENCY=LENIENT \
      $INPUT \
      OUTPUT=!{sample_and_library}.merged.bam
   '''
}

// the output of merging results in a tuple of libraries; however, each lane that is merged should be
// from the same library so just get the first value of the tuple becuase they all are the same
// downstream process depend on just a library not tuples of libraries for input
output_ch_PicardTools_MergeSamFiles_from_same_sample_and_library
   .map{ library, bams ->
      return tuple(library.get(0), bams)
   }
   .mix(output_ch_2_PicardTools_MarkDuplicates.input_ch_PicardTools_MergeSamFiles_from_same_library_PicardTools_BuildBamIndex)
	.groupTuple()
   .branch { library, bams ->
		input_ch_PicardTools_MergeSamFiles_from_same_library: bams.size() > 1

      // when grouping in values besides the first value passed, the key become wrapped in an additional tuple
      // and b/c we know that the bams are a tuple of size 1 and the downstream input requires a file 
      // that is not wrapped in a tuple, we just get the 1st element of the tuple, the file itself
		input_ch_PicardTools_BuildBamIndex: bams.size() <= 1
         return tuple(library, bams.get(0))
	}
	.set { output_ch_2_PicardTools_MarkDuplicates_PicardTools_MergeSamFiles_from_same_sample_and_library }

// merge bams from the same library with picard
process PicardTools_MergeSamFiles_from_same_library  {
   container docker_image_PicardTools

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(library), 
         file(input_bams)
      ) from output_ch_2_PicardTools_MarkDuplicates_PicardTools_MergeSamFiles_from_same_sample_and_library.input_ch_PicardTools_MergeSamFiles_from_same_library.collect()

   output:
      tuple(val(library), file("${library}.merged.bam")) into output_ch_PicardTools_MergeSamFiles_from_same_library

   shell:
   '''
   set -euo pipefail

   # add picard option prefix, 'INPUT=' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / INPUT=/g' | sed '1s/^/INPUT=/')

   java -Xmx10g -jar /picard-tools/picard.jar \
      MergeSamFiles \
      USE_THREADING=true \
      VALIDATION_STRINGENCY=LENIENT \
      $INPUT \
      OUTPUT=!{library}.merged.bam
   '''
}

// index bams with picard
process PicardTools_BuildBamIndex  {
   container docker_image_PicardTools

   publishDir params.output_dir, mode: 'copy'

   input:
      tuple(val(library), 
         file(input_bam)
      ) from output_ch_PicardTools_MergeSamFiles_from_same_library
         .mix(output_ch_2_PicardTools_MarkDuplicates_PicardTools_MergeSamFiles_from_same_sample_and_library.input_ch_PicardTools_BuildBamIndex)

   // no need for an output channel becuase this is the final stepp
   output:
      file("${library}.bam")
      file("${library}.bam.bai")

   script:
   """
   set -euo pipefail

   java -Xmx6g -jar /picard-tools/picard.jar \
      BuildBamIndex \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${library}.bam.bai

   # rename final output bam
   mv ${input_bam} ${library}.bam 
   """
}
