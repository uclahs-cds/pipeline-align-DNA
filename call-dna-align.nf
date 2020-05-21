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
      def read_group_name = "@RG" +
         "\tID:" + row.read_group_identifier + ".Seq" + row.lane +
         "\tCN:" + row.sequencing_center +
         "\tLB:" + row.library_identifier +
         "\tPL:" + row.platfrom_technology +
         "\tPU:" + row.platform_unit +
         "\tSM:" + row.sample

      return tuple(row.library_identifier,
         row.sample , 
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
      save_intermediate_files = ${params.save_intermediate_files}

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

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input: 
      tuple(val(library), 
         val(sample),
         val(lane),
         val(read_group_name), 
         path(read1_fastq),
         path(read2_fastq) 
      ) from input_samples_ch
      each file(ref) from reference
      each file(ref_dict) from reference_dict
      file(ref_idx_files) from reference_index_files.collect()

   output:
      tuple(val(library), 
         val(sample),
         val(lane),
         file("${library}-${sample}-${lane}.aligned.sam")
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
      ${library}-${sample}-${lane}.aligned.sam
   """
}

process convert_sam_to_bam {
   container "blcdsdockerregistry/samtools:1.3"

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input: 
      tuple(val(library), 
         val(sample), 
         val(lane),
         file(input_sam)
      ) from align_output_ch

   output:
      tuple (val(library), 
         val(sample), 
         val(lane), 
         file("${library}-${sample}-${lane}.aligned.bam")
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
      ${library}-${sample}-${lane}.aligned.bam
   """
}

process sort_bam  {
   container "blcdsdockerregistry/picard-tools:1.130"

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(library), 
         val(sample),
         val(lane),
         file(input_bam)
      ) from convert_sam_to_bam_output_ch
   
   output:
      tuple(val(library), 
         val(sample),
         val(lane),
         file("${library}-${sample}-${lane}.sorted.bam")
      ) into sort_bam_output_ch

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

process mark_duplicates  {
   container "blcdsdockerregistry/picard-tools:1.130"

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(library), 
         val(sample), 
         val(lane),
         file(input_bam)
      ) from sort_bam_output_ch

   output:
      tuple(val("${library}-${sample}"),
         val(library), 
         val(lane),
         file("${library}-${sample}-${lane}.mark_dup.bam")
      ) into mark_duplicates_output_ch
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

mark_duplicates_output_ch
	.groupTuple()
	.map { library_and_sample, library, lane, mark_dups_bams ->
		return tuple(library_and_sample, 
         library,
         mark_dups_bams
      )
	}
   .branch {
		to_merge_from_same_lane: it.get(2).size() > 1
		to_merge_from_same_library: it.get(2).size() <= 1
	}
	.set { mark_dups_bams_outputs }

mark_dups_bams_outputs
   .to_merge_from_same_library
   .map { library_and_sample, library, mark_dups_bams ->
      return tuple(library, mark_dups_bams)
   }
   .set { mark_dups_outputs_to_merge_bams_from_same_library }

process merge_bams_from_same_lane  {
   container "blcdsdockerregistry/picard-tools:1.130"

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(library_and_sample), 
         val(library),
         file(input_bams)
      ) from mark_dups_bams_outputs.to_merge_from_same_lane.collect()

   output:
      tuple(val(library),
         file("${library_and_sample}.merged.bam")
      ) into merge_bams_from_same_lane_output_ch

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
      OUTPUT=!{library_and_sample}.merged.bam
   '''
}

merge_bams_from_same_lane_output_ch
   .into{ merge_bams_from_same_lane_outputs_to_index_t_ch;
     merge_bams_from_same_lane_outputs_to_merge_bams_from_same_library }

merge_bams_from_same_lane_outputs_to_merge_bams_from_same_library
   .mix(mark_dups_outputs_to_merge_bams_from_same_library)
	.groupTuple()
   .branch {
		to_merge_from_same_library: it.get(1).size() > 1
		to_get_bam_index: it.get(1).size() <= 1
	}
	.set { mark_dups_and_merge_from_same_lane_outputs }

process merge_bams_from_same_library  {
   container "blcdsdockerregistry/picard-tools:1.130"

   publishDir params.output_dir, enabled: params.save_intermediate_files, mode: 'copy'

   input:
      tuple(val(library), 
         file(input_bams)
      ) from mark_dups_and_merge_from_same_lane_outputs.to_merge_from_same_library.collect()

   output:
      tuple(val(library), file("${library}.merged.bam")) into merge_bams_from_same_library_output_ch

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
      OUTPUT=!{library}.bam
   '''
}

process get_bam_index  {
   container "blcdsdockerregistry/picard-tools:1.130"

   publishDir params.output_dir, mode: 'copy'

   input:
      tuple(val(library), 
         file(input_bam)
      ) from merge_bams_from_same_library_output_ch
         .mix(mark_dups_and_merge_from_same_lane_outputs.to_get_bam_index)

   output:
      file("${input_bam.baseName.baseName}.bam.bai")
      file("${input_bam.baseName.baseName}.bam")

   script:
   """
   set -euo pipefail

   java -Xmx6g -jar /picard-tools/picard.jar \
      BuildBamIndex \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=${input_bam} \
      OUTPUT=${input_bam.baseName.baseName}.bam.bai

   # rename final output bam
   mv ${input_bam.baseName.baseName}.bam
   """
}
