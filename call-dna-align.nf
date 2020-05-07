#!/usr/bin/env nextflow 

/* 
   TODO:
      - NEEDS TESTING FOR ALIGNING MULTIPLE BAMS IN PARALLEL
      - NEEDS DOCUMENTATION 
      - NEEDS COMMENTS
*/

Channel 
   .fromFilePairs(params.paired_reads)
   .ifEmpty { error "Cannot find any reads matching: ${params.paired_reads}" }
   .set { paired_reads_ch }

Channel
   .fromPath(params.reference_fasta)
   .ifEmpty { error "Cannot find reference: ${params.reference_fasta}" }
   .set { reference }

Channel
   .fromPath(params.reference_fasta_dict)
   .ifEmpty { error "Cannot find reference dictionary: ${params.reference_fasta_dict}" }
   .set { reference_dict }

Channel
   .fromPath(params.reference_fasta_index_files)
   .ifEmpty { error "Cannot find reference index files: ${params.reference_fasta_index_files}" }
   .set { reference_index_files }

log.info """\
   =============================================
   C A L L - D N A - A L I G N  P I P E L I N E
   =============================================
   Boutros Lab

   Current Configuration:
   - input: 
      paired_reads: ${params.paired_reads}
      reference_fasta: ${params.reference_fasta}
      reference_fasta_dict: ${params.reference_fasta_dict}
      reference_fasta_index_files: ${params.reference_fasta_index_files}
      read_group_name = ${params.read_group_name}
      save_aligned_bam = ${params.save_aligned_bam}
      merge_bams = ${params.merge_bams}
   - output: 
      output_dir: ${params.output_dir}

   Tools Used:
      bwa: "blcdsdockerregistry/bwa:0.7.15"
      samtools: "blcdsdockerregistry/samtools:1.3"
      picard: "blcdsdockerregistry/picard-tools:1.130"

   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

process align {
   input: 
      tuple val(sample_name), file(paired_reads) from paired_reads_ch
      file(ref) from reference
      file(ref_dict) from reference_dict
      file(ref_idx_files) from reference_index_files.collect()

   output:
      file("${sample_name}.aligned.sam") into align_output_ch

   container = "blcdsdockerregistry/bwa:0.7.15"

   script:
   """
   set -euo pipefail

   bwa \
      mem \
      -t 32 \
      -M \
      -R "${params.read_group_name}" \
      ${ref} \
      ${paired_reads[0]} \
      ${paired_reads[1]} > \
      ${sample_name}.aligned.sam
   """
}

process convert_sam_to_bam {
   container "blcdsdockerregistry/samtools:1.3"

   publishDir params.output_dir, enabled: params.save_aligned_bam

   input: 
      file(input_sam) from align_output_ch

   output:
      file("${input_sam.baseName}.bam") into convert_sam_to_bam_output_ch

   script:
   """
   set -euo pipefail

   samtools \
      view \
      -@ 4 \
      -S \
      -b \
      ${input_sam} > \
      ${input_sam.baseName}.bam
   """
}

process sort_bam  {
   container "blcdsdockerregistry/picard-tools:1.130"

   input:
      file(input_bam) from convert_sam_to_bam_output_ch
   
   output:
      file("${input_bam.baseName}.sorted.bam") into sort_bam_output_ch

   script:
   """
   set -euo pipefail

   java -Xmx14g -jar /picard-tools/picard.jar \
      SortSam \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${input_bam.baseName}.sorted.bam \
      SORT_ORDER=coordinate
   """
}

process mark_duplicates  {
   container "blcdsdockerregistry/picard-tools:1.130"

   input:
      file(input_bam) from sort_bam_output_ch

   output:
      file("${input_bam.baseName}.mark_dup.bam") into mark_duplicates_output_ch
      file("${input_bam.baseName}.mark_dup.bam.metrics") into mark_duplicates_metrics_ch

   script:
   """
   set -euo pipefail

   java -Xmx10g -jar /picard-tools/picard.jar \
      MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${input_bam.baseName}.mark_dup.bam \
      METRICS_FILE=${input_bam.baseName}.mark_dup.bam.metrics \
      PROGRAM_RECORD_ID=MarkDuplicates
   """
}

// get the number of aligned bams and 
(merge_bams_input_ch, get_bam_index_input_ch) = ( params.merge_bams
   ? [  mark_duplicates_output_ch, Channel.empty() ]
   : [  Channel.empty(), mark_duplicates_output_ch ] )

process merge_bams  {
   container "blcdsdockerregistry/picard-tools:1.130"

   input:
      file(input_bams) from merge_bams_input_ch.collect()

   output:
      file("merged.bam") into merge_bams_output_ch

   shell:
   '''
   set -euo pipefail

   # add picard option prefix, 'INPUT=' to each input bam
   declare -r INPUT=$(echo '!{input_bams}' | sed -e 's/ / INPUT=/g' | sed '1s/^/INPUT=/')

   java -Xmx6g -jar /picard-tools/picard.jar \
      MergeSamFiles \
      USE_THREADING=true \
      VALIDATION_STRINGENCY=LENIENT \
      $INPUT \
      OUTPUT=merged.bam
   '''
}

process get_bam_index  {
   container "blcdsdockerregistry/picard-tools:1.130"

   publishDir params.output_dir

   input:
      file(input_bam) from get_bam_index_input_ch.mix(merge_bams_output_ch)

   output:
      file("${input_bam.baseName}.bam.bai") into final_bam_index
      file("${input_bam.baseName}.bam") into final_bam

   script:
   """
   set -euo pipefail

   java -Xmx6g -jar /picard-tools/picard.jar \
      BuildBamIndex \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${input_bam.baseName}.bam.bai
   """
}
