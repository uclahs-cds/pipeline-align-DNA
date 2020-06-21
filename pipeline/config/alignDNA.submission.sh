#!/bin/bash -l

# generate SBATCH/Slurm options
##########################################################################################
#SBATCH --exclusive # reserve an entire node
#SBATCH --partition=midmem # worker node partition type (e.g: midmem or execute)
#SBATCH -J alignDNA-a-partial- # name of the run in the queue 
#SBATCH -e alignDNA-a-partial-.error # name of the stderror file 
#SBATCH -o alignDNA-a-partial-.log # name of stdout file 
#SBATCH --mail-user=bcarlin # Where to send email
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
    declare -r pipeline_run_name=$3
    declare -r email=$4

    # get the nextflow run command and execute the pipeline
    execute_nextflow_run_command $nextflow_script $config_file
}

# 1st input ($1): nextflow script (.nf)
# 2nd  input ($2): config file (.config)
# 3rd input ($3): pipeline run name (alignDNA-<dataset name>)
# 4th input ($4): email (<AD>@mednet.ucla.edu) 
main $1 $2 $3 $4