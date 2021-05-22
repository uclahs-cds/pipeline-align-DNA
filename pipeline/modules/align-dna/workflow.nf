include { validate_file } from '../validation.nf'
include { align_DNA_BWA_MEM2; align_DNA_HISAT2 } from './fastq2bam.nf'
include { PicardTools_SortSam; PicardTools_MarkDuplicates; PicardTools_BuildBamIndex } from './picardtools.nf'
include { Generate_Sha512sum } from './checksum.nf'

workflow aligndna {
   take:
      ich_samples
      ich_reference_fasta
      ich_reference_index_files
   main:
      if (params.aligner == "BWA-MEM2") {
         align_DNA_BWA_MEM2(
            ich_samples,
            ich_reference_fasta,
            ich_reference_index_files.collect()
            )
         PicardTools_SortSam(align_DNA_BWA_MEM2.out.bam)
      }	 
      else if (params.aligner == "HISAT2") { 
         align_DNA_HISAT2(
            ich_samples,
            ich_reference_fasta,
            ich_reference_index_files.collect()
            )
	 PicardTools_SortSam(align_DNA_HISAT2.out.bam)	 
      } 
      PicardTools_MarkDuplicates(PicardTools_SortSam.out.bam.collect())
      PicardTools_BuildBamIndex(PicardTools_MarkDuplicates.out.bam)
      Generate_Sha512sum(PicardTools_BuildBamIndex.out.bai.mix(PicardTools_MarkDuplicates.out.bam))
      validate_file(
         PicardTools_MarkDuplicates.out.bam.mix(
            PicardTools_BuildBamIndex.out.bai,
            Channel.from(params.temp_dir, params.output_dir)
            )
         )
}
