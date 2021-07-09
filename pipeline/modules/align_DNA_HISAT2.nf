// The follwing process runs both alignment and SAM conversion in the same process.
// While this is normally considered to go against the best practices for processes,
// here it actually saves cost, time, and memory to directly pipe the output into 
// samtools due to the large size of the uncompressed SAM files.

include { validate_file; validate_file as validate_output_file } from './run_validate.nf'
include { PicardTools_SortSam } from './sort_bam_picardtools.nf'
include { PicardTools_MarkDuplicates } from './mark_duplicate_picardtools.nf'
include { PicardTools_BuildBamIndex } from './index_bam_picardtools.nf'
include { Generate_Sha512sum } from './check_512sum.nf'

process align_DNA_HISAT2 {
   container params.docker_image_hisat2_and_samtools
   publishDir path: params.bam_output_dir,
      enabled: params.save_intermediate_files,
      pattern: "*.bam",
      mode: 'copy'

   publishDir path: params.log_output_dir,
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "align_DNA_HISAT2/${file(read1_fastq).getSimpleName()}/${library}-${lane}.log${file(it).getName()}" }

   // use "each" so the the reference files are passed through for each fastq pair alignment 
   input: 
      tuple(val(library),
         val(header), 
         val(lane), 
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

   hisat2 \
      -p ${task.cpus} \
      --rg-id ${header.read_group_identifier}.Seq${header.lane} \
      --rg LB:${header.library_identifier} \
      --rg PL:${header.platform_technology} \
      --rg PU:${header.platform_unit} \
      --rg SM:${header.sample} \
      --no-spliced-alignment \
      -x ${file(ich_reference_index_files.get(0).baseName).baseName} \
      -1 ${read1_fastq} \
      -2 ${read2_fastq} | \
   samtools \
      view \
      -@ ${task.cpus} \
      -S \
      -b > \
      ${library}-${lane}.bam
   """
   }

workflow align_DNA_HISAT2_workflow {
   aligner_output_dir = "${params.bam_output_dir}/HISAT2-2.2.1"
   take:
      ich_samples
      ich_samples_validate
      ich_reference_fasta
      ich_reference_index_files
   main:
      validate_file(ich_samples_validate.mix(
         ich_reference_fasta,
         ich_reference_index_files
         ))
      align_DNA_HISAT2(
         ich_samples,
         ich_reference_fasta,
         ich_reference_index_files.collect()
         )
      PicardTools_SortSam(align_DNA_HISAT2.out.bam)
      PicardTools_MarkDuplicates(PicardTools_SortSam.out.bam.collect(), aligner_output_dir)
      PicardTools_BuildBamIndex(PicardTools_MarkDuplicates.out.bam, aligner_output_dir)
      Generate_Sha512sum(PicardTools_BuildBamIndex.out.bai.mix(PicardTools_MarkDuplicates.out.bam), aligner_output_dir)
      validate_output_file(
         PicardTools_MarkDuplicates.out.bam.mix(
            PicardTools_BuildBamIndex.out.bai,
            Channel.from(params.temp_dir, params.output_dir)
            )
         )
   }

