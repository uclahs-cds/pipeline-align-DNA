process run_validate {
   container params.docker_image_validate_params

   publishDir path: "${log_output_dir}/${task.process.split(':')[1].replace('_', '-')}/${task.index}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "log${file(it).getName()}" }

   input:
   path(file_to_validate)
   val(log_output_dir)

   output:
   path(".command.*")
   path("input_validation.txt"), emit: val_file

   script:
   """
   set -euo pipefail

   python -m validate -t file-input ${file_to_validate} > "input_validation.txt"
   """
   }
