#!/bin/bash -l

# generate SBATCH/Slurm options
function generate_sbatch_options() {
    # get the pipeline run name and email
    declare -r pipeline_run_name=$1
    declare -r email=$2

    declare -r sbatch_options="
#SBATCH --exclusive # reserve an entire node
#SBATCH --partition=midmem # worker node partition type (e.g: midmem or execute)
#SBATCH -J ${pipeline_run_name} # name of the run in the queue 
#SBATCH -e ${pipeline_run_name}.error # name of the stderror file 
#SBATCH -o ${pipeline_run_name}.log # name of stdout file 
#SBATCH --mail-user=${email} # Where to send email
#SBATCH --mail-type=ALL  # Mail events (NONE, BEGIN, END, FAIL, ALL)
"

    # return the options
    echo $sbatch_options
}

# create temp directory for /work and intermediates file
function create_work_dir() {
    # get the pipeline run name 
    declare -r pipeline_run_name=$1

    # create temp directory for /work and intermediates file 
    declare -a work_dir=$(mktemp -d /scratch/$pipeline_run_name-temp-XXXXXX)

    # return path of the work dir
    echo $work_dir
}

# create the nextflow run command
function execute_nextflow_run_command() {
    # get the nextflow script, configuration file and work dir
    declare -r nextflow_script=$1
    declare -r config_file=$2
    declare -r work_dir=$3

    # nextflow run command
    declare -r nextflow_run_command="nextflow run $nextflow_script -config $config_file -work-dir $work_dir"

    # execute nextflow run command
    echo $nextflow_run_command
    eval $nextflow_run_command
}

function main() {
    # get the pipeline run name, nextflow script, configuration file and email
    declare -r nextflow_script=$1
    declare -r config_file=$2
    declare -r pipeline_run_name=$3
    declare -r email=$4

    # generate SBATCH/Slurm options
    generate_sbatch_options $pipeline_run_name $email

    # create temp directory for /work and intermediates file 
    declare -a work_dir=$(create_work_dir $pipeline_run_name)

    # get the nextflow run command and execute the pipeline
    execute_nextflow_run_command $nextflow_script $config_file $work_dir

    # clean the work directory
    rm -r $work_dir
}

# 1st input ($1): nextflow script (.nf)
# 2nd  input ($2): config file (.config)
# 3rd input ($3): pipeline run name (alignDNA-<dataset name>)
# 4th input ($4): email (<AD>@mednet.ucla.edu) 
main $1 $2 $3 $4