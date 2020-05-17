Channel
   .fromPath(params.input_csv)
   .ifEmpty { error "Cannot find input csv: ${params.input_csv}" }
   .splitCsv(header:true)
   .map { row -> 
      return tuple(row.bam_1, row.bam_2)
   }
   .set { input_bams_ch }
   //.subscribe { println it }

log.info """\
   ==========================================
   C O M P A R E - B A M S - P I P E L I N E
   ==========================================
   Boutros Lab

   Current Configuration:
   - input: 
      input_csv: ${params.input_csv}
   - output: 
      output_dir: ${params.output_dir}
   - options:
      get_differences: ${params.get_differences}
   
   Tools Used:
      jvarkit's cmpbams: "blcdsdockerregistry/jvarkit-cmpbams:1.0"
    
   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

process compare_bams {
   container "blcdsdockerregistry/jvarkit-cmpbams:1.0"

   publishDir params.output_dir, mode: 'move'

   input: 
      tuple(path(bam_1), path(bam_2)) from input_bams_ch

   output:
     file("${bam_2.baseName}-${bam_2.baseName}-cmp.txt") into compare_bams_output_ch

   script:
   """
   set -euo pipefail

   java -Xmx24g -jar /jvarkit-cmpbams/cmpbams.jar \
      ${bam_1} \
      ${bam_2} \
      > ${bam_2.baseName}-${bam_2.baseName}-cmp.txt
   """
}

process find_different_reads {
   container "ubuntu:latest"

   publishDir params.output_dir, mode: 'move' 

   input: 
      file(bam_1_bam_2_comparison) from compare_bams_output_ch

   output:
     file("${bam_2.baseName}-${bam_2.baseName}-cmp.txt") into find_different_reads_output_ch

   when: 
      params.get_differences == true

   script:
   """
   cat ${bam_1_bam_2_comparison} \
      | grep -v "EQ" \
      >  ${bam_1_bam_2_comparison}.diff.txt
   """
}
