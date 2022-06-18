#!/bin/bash 
#SBATCH -J test-memory-samtools-sort
#SBATCH -o test-memory-samtools-sort.out
#SBATCH -e test-memory-samtools-sort.error
#SBATCH --partition="F16"
#SBATCH --exclusive
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=jarbet@mednet.ucla.edu
#SBATCH --profile=Task
#SBATCH --acctg-freq=1

task_cpus=12;
path_out=/hot/software/pipeline/pipeline-align-DNA/Nextflow/development/unreleased/jarbet-memory-samtools-sort;
path_input_bam=/hot/software/pipeline/pipeline-align-DNA/Nextflow/development/unreleased/jarbet-undef-bam-output-dir/BWA-MEM2/a_mini_n1/align-DNA-8.0.0/a_mini_n1/BWA-MEM2-2.2.1/output;
input_bam=a_mini_n1.bam;
sort_order="-n";
bam_output_filename=${input_bam}.sorted;

docker run -u $(id -u):$(id -g) \
    -v ${path_input_bam}:${path_input_bam} \
    -v ${path_out}:${path_out} \
    -w ${path_out} \
    blcdsdockerregistry/samtools:1.15.1 \
    samtools sort \
        -@ ${task_cpus} \
        -O bam \
        -o ${path_out}/${bam_output_filename} \
        ${sort_order} \
        ${path_input_bam}/${input_bam};

