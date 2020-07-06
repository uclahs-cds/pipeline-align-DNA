# pipeline-align-DNA

For pipeline documentation, please refer [here](https://uclahs.box.com/s/kl4pacq332bprpe9lnfams0l8vglmg30)

For docker images, please view our Dockerhub repository [here](https://hub.docker.com/orgs/blcdsdockerregistry/repositories)

## How to run the pipeline
**The pipeline should be run WITH A SINGLE SAMPLE AT TIME. Otherwise resource allocation and Nextflow errors could cause the pipeline to fail**

On your submitter node run:
sbatch or qsub /path/to/align-DNA.submission.sh /path/to/align-DNA.nf /path/to/align-DNA.config
