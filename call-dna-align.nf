#!/usr/bin/env nextflow 

/* 
   TODO:
      - NEEDS TESTING FOR ALIGNING MULTIPLE BAMS IN PARALLEL
      - NEEDS DOCUMENTATION 
      - NEEDS COMMENTS
*/

Channel
   .fromPath(params.input_csv)
   .ifEmpty { error "Cannot find input csv: ${params.input_csv}" }
   .splitCsv(header:true)
   .map { row -> 
      def sample_name = row.read_group_identifier + "-" +
         row.sequencing_center + "-" +
         row.library_identifier + "-" +
         row.platfrom_technology + "-" +
         row.platform_unit + "-" +
         row.sample 

      def read_group_name = "@RG" +
         "\tID:" + row.read_group_identifier + ".Seq" + row.lane +
         "\tCN:" + row.sequencing_center +
         "\tLB:" + row.library_identifier +
         "\tPL:" + row.platfrom_technology +
         "\tPU:" + row.platform_unit +
         "\tSM:" + row.sample

      return tuple(sample_name, 
         row.lane,
         read_group_name,
         row.read1_fastq,
         row.read2_fastq
      )
   }
   .set { input_samples_ch }


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
      input_csv: ${params.input_csv}
      reference_fasta: ${params.reference_fasta}
      reference_fasta_dict: ${params.reference_fasta_dict}
      reference_fasta_index_files: ${params.reference_fasta_index_files}
   - output: 
      output_dir: ${params.output_dir}
   - options:
      save_aligned_bam = ${params.save_aligned_bam}

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
   container "blcdsdockerregistry/bwa:0.7.15"

   input: 
      tuple (val(sample_name), 
         val(lane),
         val(read_group_name), 
         path(read1_fastq),
         path(read2_fastq) 
      ) from input_samples_ch
      each file(ref) from reference
      each file(ref_dict) from reference_dict
      file(ref_idx_files) from reference_index_files.collect()

   output:
      tuple(val(sample_name), 
         val(lane),
         file("${sample_name}.lane-${lane}.aligned.sam")
      ) into align_output_ch

   script:
   """
   set -euo pipefail

   bwa \
      mem \
      -t 32 \
      -M \
      -R "${read_group_name}" \
      ${ref} \
      ${read1_fastq} \
      ${read2_fastq} > \
      ${sample_name}.lane-${lane}.aligned.sam
   """
}

process convert_sam_to_bam {
   container "blcdsdockerregistry/samtools:1.3"

   publishDir params.output_dir, enabled: params.save_aligned_bam

   input: 
      tuple(val(sample_name), 
         val(lane),
         file(input_sam)
      ) from align_output_ch

   output:
      tuple(val(sample_name), 
         val(lane), 
         file("${sample_name}.lane-${lane}.converted.bam")
      ) into convert_sam_to_bam_output_ch

   script:
   """
   set -euo pipefail

   samtools \
      view \
      -@ 4 \
      -S \
      -b \
      ${input_sam} > \
      ${sample_name}.lane-${lane}.converted.bam
   """
}

process sort_bam  {
   container "blcdsdockerregistry/picard-tools:1.130"

   input:
      tuple(val(sample_name), 
      val(lane),
      file(input_bam)
   ) from convert_sam_to_bam_output_ch
   
   output:
      tuple(val(sample_name),
         val(lane),
         file("${sample_name}.lane-${lane}.sorted.bam")
      ) into sort_bam_output_ch

   script:
   """
   set -euo pipefail

   java -Xmx14g -jar /picard-tools/picard.jar \
      SortSam \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${sample_name}.lane-${lane}.sorted.bam \
      SORT_ORDER=coordinate
   """
}

process mark_duplicates  {
   container "blcdsdockerregistry/picard-tools:1.130"

   input:
      tuple(val(sample_name), 
         val(lane),
         file(input_bam)
      ) from sort_bam_output_ch

   output:
      tuple(
         val(sample_name),
         val(lane),
         file("${sample_name}.lane-${lane}.mark_dup.bam")
      ) into mark_duplicates_output_ch
      file("${sample_name}.lane-${lane}.mark_dup.bam.metrics") into mark_duplicates_metrics_ch

   script:
   """
   set -euo pipefail

   java -Xmx10g -jar /picard-tools/picard.jar \
      MarkDuplicates \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${sample_name}.lane-${lane}.mark_dup.bam \
      METRICS_FILE=${sample_name}.lane-${lane}.mark_dup.bam.metrics \
      PROGRAM_RECORD_ID=MarkDuplicates
   """
}

mark_duplicates_output_ch
	.groupTuple()
	.map { sn, lane, mark_dups_bams ->
		def sample_name = sn
		return tuple(sample_name, mark_dups_bams, mark_dups_bams.size())
	}
   .branch {
		merge: it.get(2) > 1
		index: it.get(2) <= 1
	}
	.set { mark_dups_bams_outputs }

mark_dups_bams_outputs.merge.set { merge_bams_input_ch }
mark_dups_bams_outputs
   .index
   .map { sample_name, mark_dups_bams, size ->
      return tuple(sample_name, mark_dups_bams)
   }
   .set { get_bam_index_input_ch }


process merge_bams  {
   container "blcdsdockerregistry/picard-tools:1.130"

   input:
      tuple(val(sample_name), file(input_bams), val(input_bams_size)) from merge_bams_input_ch

   output:
      tuple(val(sample_name), file("${sample_name}.merged.bam")) into merge_bams_output_ch

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
      OUTPUT=!{sample_name}.merged.bam
   '''
}

process get_bam_index  {
   container "blcdsdockerregistry/picard-tools:1.130"

   publishDir params.output_dir

   input:
      tuple(
         val(sample_name), 
         file(input_bam)
      ) from get_bam_index_input_ch.mix(merge_bams_output_ch)

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
