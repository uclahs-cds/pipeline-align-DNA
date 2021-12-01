# Changelog

All notable changes to this project will be documented in this file.


The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- Changed names of midmem.config and execute.config into F72.config and M64.config respectively.
- Rename bug report to "Issue report" and remove old node names from it

## [7.2.1] - 2021-10-28
### Changed
- Fix potential Spark temp directory permissions issue

## [7.2.0] - 2021-10-01
### Added
- GPL2 License added
- MarkDuplicatesSpark process added as an option
### Changed
- Removed explicit index creation process and enabled option for MarkDuplicate process to create index
- Allow CPU and memory allocation to dictate parallelization rather than maxForks

## [7.1.0] - 2021-07-29
### Added
- HISAT2 aligner functionality and the option to run either BWA-MEM2/HISAT2 or both at once. The default aligner is BWA-MEM2.
- A python script to generate config files from command line.
### Changed
- Update config file to process inputs for each aligner separately. Old config files still work and BWA-MEM2 will be run as usual.
- #112 Update BWA-MEM2 and SAMtools docker to SAMtools 1.12
- #121 Update version information in the main script
- #126 Update output directory structure
- Process names standardized
- #128 Use explicit tab delimiters to ensure proper program tagging
- Updated validation docker image to v2.1.5
### Removed
- Dockerfiles for BWA-MEM2, jvarkit-cmpbams, and Picard removed and moved to their own separate repositories (docker-BWA-MEM2, docker-jvarkit-cmpbams, and docker-Picard, respectively).

## [7.0.3] - 2021-04-15
### Changed
- #61 Update validation to 2.1.0.
- #76 Update version documentation and manifest.
- #78 #81 update resources setting for alignment, sort, and markduplicate
- #79 Update CHANGELOG.md to reflect Keep a Changelog format.
- #82 Save outputs in directories based on FASTQ library/lane #2.
- #83 Rename main workflow module.
- #88 node specific configs are not included properly
- #89 docker permission is not set properly
- #90 Fixed dockerfiles to pass dockerfilelint


## [7.0.2] - 2021-03-15
### Changed
- #70 Fixes crash related to checking default node CPU and memory configurations
- #67 Run docker with group permissions of the user executing the pipeline
- #65 Check for write permission on output directories before executing


## [7.0.1] - 2021-03-11 [YANKED]
### Changed
- #67 Run docker with group permissions of the user executing the pipeline
- #65 Check for write permission on output directories before executing


## [7.0.0] - 2021-03-02
### Changed
- #43 Port to DSL2


## [6.1.0] - 2021-01-18
### Added
- A small pipeline to generate the reference genome index files. This is a separate nextflow script from the main pipeline script.


## [6.0.2] - 2021-01-13
### Fixed
- Processes in a Docker container are executed as the user automatically instead of root.


## [6.0.1] - 2020-12-08
### Changed
- Process name for alignment is changed to align\_BWA\_mem\_convert\_SAM\_to\_BAM\_samtools to be readable
- Gave sudo to Docker when running on the sge cluster

### Fixed
- #31 Error: Unknown method invocation `includeConfig`
- #32 Error: Please specify the disease\_id, patient\_id, dataset\_id, sample\_id, analyte, and technology in the config file


## [6.0.0] - 2020-11-30
### Added
- Validation scripts are fully implemented. The pipeline will stop if invalid input/output files detected, e.g., files not found are wrong file type
- Enabled input and output directly from and to the Boutros Lab data storage

### Changed
- Simplified the config file opened to users with only essensial parameters included


## [5.0.0] - 2020-11-09
### Changed
- bwa-mem2 is upgraded to v2.1. It provides a smaller indexed genome and lower cpu usage comparing to the previous version v2.0.


## [4.0.0-beta] - 2020-10-15
### Added
- Nextflowization of align-DNA pipeline
- Dynamic resource allocation

### Changed
- Version tool updates (BWA 0.7.17, Samtools 1.10, Picard Tools 2.23.3)


## [0.0.1] - 2020-10-08
### Added
- Initial Release
