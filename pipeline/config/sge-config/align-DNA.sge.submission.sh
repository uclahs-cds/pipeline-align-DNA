#!/bin/bash -l

# generate SGE options
#########################################################################################
#$ -S /bin/bash 
#$ -l slot_type=<partition_type>,h_vmem=<amount_of_mem>,h=<ip-for_specific_worker_node>
#$ -cwd # place the log files/output in the current working directory 
#$ -N <pipeline_run_name> # name of the run in the queue 
#$ -e <pipeline_run_name> # name of the stderror file 
#$ -o <pipeline_run_name> # name of stdout file 
#$ -m abe # circumstances upon which to send an email: abort, beginning, end
#$ -pe smpslots <partition slots> # parameter for the cloud usage per P Tung 
#$ -M <email address> # e-mail address where messages are sent 
#$ -r no # do not rerun the job if it fails 
##########################################################################################

# create the nextflow run command
function execute_nextflow_run_command() {
    # get the nextflow script, configuration file and work dir
    declare -r nextflow_script=$1
    declare -r config_file=$2

    # nextflow run command
    declare -r nextflow_run_command="nextflow run $nextflow_script -config $config_file"

    # output and then execute nextflow run command
    eval $nextflow_run_command
}

function main() {
    # get the pipeline run name, nextflow script, configuration file and email
    declare -r nextflow_script=$1
    declare -r config_file=$2

    # get the nextflow run command and execute the pipeline
    execute_nextflow_run_command $nextflow_script $config_file
}

# 1st input ($1): nextflow script (.nf)
# 2nd  input ($2): config file (.config)
main $1 $2