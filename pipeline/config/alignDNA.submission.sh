#!/bin/bash -l

# generate SBATCH/Slurm options
##########################################################################################
#SBATCH --exclusive # reserve an entire node
#SBATCH --partition=midmem # worker node partition type (e.g: midmem or execute)
#SBATCH -J <pipeline_run_name> # name of the run in the queue 
#SBATCH -e <pipeline_run_name>.error # name of the stderror file 
#SBATCH -o <pipeline_run_name>.log # name of stdout file 
#SBATCH --mail-user=<email address> # Where to send email
#SBATCH --mail-type=ALL  # Mail events (NONE, BEGIN, END, FAIL, ALL)
##########################################################################################

# create the nextflow run command
function execute_nextflow_run_command() {
    # get the nextflow script, configuration file and work dir
    declare -r nextflow_script=$1
    declare -r config_file=$2

    # nextflow run command
    declare -r nextflow_run_command="nextflow run $nextflow_script -config $config_file"

    # output and then execute nextflow run command
    echo $nextflow_run_command
    eval $nextflow_run_command
}

function main() {
    # make sure you are on local disk
    cd /scratch

    # get the pipeline run name, nextflow script, configuration file and email
    declare -r nextflow_script=$1
    declare -r config_file=$2

    # get the nextflow run command and execute the pipeline
    execute_nextflow_run_command $nextflow_script $config_file
}

# 1st input ($1): nextflow script (.nf)
# 2nd  input ($2): config file (.config)
main $1 $2