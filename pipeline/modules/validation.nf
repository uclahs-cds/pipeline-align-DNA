process validate_file {
   container params.docker_image_validate_params

   input:
   path file_to_validate

   script:
   """
   set -euo pipefail

   python -m validate -t file-input ${file_to_validate} 1> /dev/null
   """
   }
