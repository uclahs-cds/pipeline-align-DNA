// produce checksum for the bam and bam index
process Generate_Sha512sum {    
   container params.docker_image_sha512sum

   publishDir path: params.bam_output_dir, mode: 'copy'

   input:
      file(input_file)

   output:
      file("${input_file.getName()}.sha512")

   script:
   """
   set -euo pipefail

   sha512sum ${input_file} > ${input_file.getName()}.sha512
   """
   }
