include { validate_file } from '../validation.nf'
include { Align_BWA_mem_convert_SAM_to_BAM_samtools } from './samtools.nf'
include { PicardTools_SortSam; PicardTools_MarkDuplicates; PicardTools_BuildBamIndex } from './picardtools.nf'
include { Generate_Sha512sum } from './checksum.nf'

workflow aligndna {
   take:
      ich_samples
      ich_reference_fasta
      ich_reference_index_files
   main:
      Align_BWA_mem_convert_SAM_to_BAM_samtools(
         ich_samples,
         ich_reference_fasta,
         ich_reference_index_files.collect()
      )
      PicardTools_SortSam(Align_BWA_mem_convert_SAM_to_BAM_samtools.out)
      PicardTools_MarkDuplicates(PicardTools_SortSam.out.collect())
      PicardTools_BuildBamIndex(PicardTools_MarkDuplicates.out)
      Generate_Sha512sum(PicardTools_BuildBamIndex.out.mix(PicardTools_MarkDuplicates.out))
      validate_file(
         PicardTools_MarkDuplicates.out.mix(
            Generate_Sha512sum.out,
            PicardTools_BuildBamIndex.out,
            Channel.from(params.temp_dir, params.output_dir)
         )
      )
}