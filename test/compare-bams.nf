def jvarkit_docker_image = "blcdsdockerregistry/jvarkit-cmpbams:1.0"

Channel
   .fromPath(params.input_csv)
   .ifEmpty { error "Cannot find input csv: ${params.input_csv}" }
   .splitCsv(header:true)
   .map { row -> 
      return tuple(row.bam_1, row.bam_2)
   }
   .set { input_ch_input_bams }

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
   
   Tools Used:
      jvarkit's cmpbams: ${jvarkit_docker_image}
    
   ------------------------------------
   Starting workflow...
   ------------------------------------
   """
   .stripIndent()

process jvarkit_compare_bams {
   container jvarkit_docker_image

   publishDir params.output_dir, mode: 'copy'

   input: 
      tuple(path(bam_1), path(bam_2)) from input_ch_input_bams

   output:
     file("${bam_1.baseName}-${bam_2.baseName}-diff.txt")

   script:
   """
   set -euo pipefail

   java -Xmx24g -jar /jvarkit-cmpbams/cmpbams.jar \
      ${bam_1} \
      ${bam_2} \
      | grep -v "EQ" \
      > ${bam_1.baseName}-${bam_2.baseName}-diff.txt
   """
}
