# vortex_knight

#### Description

`vortex_knight` is a taxonomic profiling nextflow workflow, utilising an ensemble of taxonomic profilers. 

---
# Requirements

The easiest way to handle dependencies is via Singularity/Docker containers. Alternatively, conda environments, software module systems or native installations can be used.

## Preprocessing

Preprocessing and QA is done with `bbmap`, `fastqc`, and `multiqc`.


---
# Usage
## Cloud-based Workflow Manager (CloWM)
This workflow will be available on the `CloWM` platform (coming soon).

## Command-Line Interface (CLI)
The workflow run is controlled by environment-specific parameters (see [run.config](https://github.com/cschu/vortex_knight/blob/main/config/run.config)) and study-specific parameters (see [params.yml](https://github.com/cschu/vortex_knight/blob/main/config/params.yml)). The parameters in the `params.yml` can be specified on the command line as well.

You can either clone this repository from GitHub and run it as follows
```
git clone https://github.com/cschu/vortex_knight.git
nextflow run /path/to/vortex_knight [-resume] -c /path/to/run.config -params-file /path/to/params.yml
```

Or, you can have nextflow pull it from github and run it from the `$HOME/.nextflow` directory.
```
nextflow run cschu/vortex_knight [-resume] -c /path/to/run.config -params-file /path/to/params.yml
```

## Input files
Fastq files are supported and can be either uncompressed (but shouldn't be!) or compressed with `gzip` or `bzip2`. Sample data can be arranged in one single input directory ("flat") or as one directory per sample ("tree").

Mates 1 and 2 can be specified with suffixes `_[12]`, `_R[12]`, `.[12]`, `.R[12]`. Lane IDs or other read id modifiers have to precede the mate identifier. Files with names not containing either of those patterns will be assigned to be single-ended. Samples consisting of both single and paired end files are assumed to be paired end with all single end files being orphans (quality control survivors). 


### All files in one directory -- "flat"
Files in the input directory must have perfectly matching prefixes in order to be associated as belonging to the same sample. Orphans belonging to the same sample as a paired-end file pair must contain the same sample prefix as well as the string `.singles` (preceding any `R1`/`R2` suffix and the fastq suffix.)

### Per-sample input subdirectories -- "tree"
All files in a sample directory will be associated with the name of the sample folder. Paired-end mate files need to have matching prefixes. 
