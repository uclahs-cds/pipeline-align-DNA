DON'T USE
workflow aligndna {
    main:
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
    emit:
      ich_samples
}