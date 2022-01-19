process run_validate {
   container params.docker_image_validate_params

   publishDir path: "${aligner_validation_dir}",
      pattern: ".command.*",
      mode: "copy",
      saveAs: { "${task.process}/${task.process}-${task.index}/log${file(it).getName()}" }

   input:
   path(file_to_validate)
   val(aligner_validation_dir)

   output:
   path(".command.*")
   path("input_validation.txt"), emit: val_file

   script:
   """
   set -euo pipefail

   python -m validate -t file-input ${file_to_validate} > 'input_validation.txt'
   """
   }
