docker_image_validate_params = "blcdsdockerregistry/validate:1.0.0"

log.info """\
------------------------------------
         V A L I D A T I O N
-----------------------------------
"""

process validate_file {
   //container docker_image_validate_params

   input:
   path file_to_validate

   script:
   """
   set -euo pipefail

   #python -m validate -t file-input ${file_to_validate}
   echo 'hello'
   """
}