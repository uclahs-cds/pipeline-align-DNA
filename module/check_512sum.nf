// produce checksum for the bam and bam index
process generate_sha512sum {    
   container params.docker_image_validate

   publishDir path: "${checksum_output_dir}",
      mode: "copy",
      pattern: "*.sha512",
      saveAs: { filename -> (filename.endsWith(".bai.sha512") && !filename.endsWith(".bam.bai.sha512")) ? "${file(file(filename).baseName).baseName}.bam.bai.sha512" : "${filename}"}

   input:
      file(input_file)
      val(checksum_output_dir)

   output:
      file("${input_file.getName()}.sha512")

   script:
   """
   set -euo pipefail

   sha512sum ${input_file} > ${input_file.getName()}.sha512
   """
   }
