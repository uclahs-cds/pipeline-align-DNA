include { validate_file } from './validation.nf'

workflow validate_ichsamples {
   main:
      Channel
         .fromPath(params.reference_fasta)
         .ifEmpty { error "Cannot find reference: ${params.reference_fasta}" }
         .set { ich_reference_fasta }

      Channel
         .fromPath(params.reference_fasta_index_files)
         .ifEmpty { error "Cannot find reference index files: ${params.reference_fasta_index_files}" }
         .set { ich_reference_index_files }

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
         .set{ ich_samples }
      
      ich_samples
         .flatMap { library, lane, read_group_name, read1_fastq, read2_fastq ->
            [read1_fastq, read2_fastq]
         }
         .set { ich_samples_validate }

      validate_file(ich_samples_validate.mix(
         ich_reference_fasta,
         ich_reference_index_files
      ))

   emit:
      ich_samples = ich_samples
      // TODO this is to expose library, lane, etc.
      // there's gotta be a better way ...
      //ich_samples_validate = ich_samples_validate
      ich_reference_fasta = ich_reference_fasta
      ich_reference_index_files = ich_reference_index_files//.collect()
}