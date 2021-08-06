# vortex_knight

### Installation

1. Clone the repo from GitHub.

```
git clone https://github.com/cschu/vortex_knight.git
```

2. Create a conda environment with NextFlow, e.g. by using the provided `environment.yml`.

```
cd vortex_knight
conda env create -f environment.yml
conda activate vortex_knight
```

3. Make a copy of the `nextflow/run.config` file and adjust it to your environment.

4. Run the pipeline 

``` 
nextflow run /path/to/vortex_knight/nextflow/vknight.nf --input_dir /path/to/input_files --bam_input (true|false) --output_dir /path/to/output_dir -c /path/to/run.config"
```

* `input_dir` should be a folder with bam files or with gzipped fastq files. For fastq files, individual samples should be separated into individual folders.
* `bam_input` must be `true` if the input are bam files and `false` if not.
* `output_dir` is `vknight_out` in the local directory by default.
* For submission to a cluster, nextflow itself requires at least `5GB` of memory.

5. Outputs

The output folder contains:

* one subdirectory `otu_tables` containing the summarised `mapseq` otu tables
* a subdirectory per sample (named `<sample>`) with
  * the kraken2 report `<sample>.kraken2_report.txt`
  * the library size `<sample>.libsize.txt`
  * the mOTUs report `<sample>.motus.txt`
  * pathseq output
    - `<sample>.pathseq.bam`
    - `<sample>.pathseq.bam.sgi`
    - `<sample>.pathseq.score_metrics`
    - `<sample>.pathseq.scores`

**Note that by default, all files in the output folder are symlinks into the work dir! Before you delete the work dir, ensure you have dereferenced copies. Alternatively, change the publishMode to `copy` or `link` (if the target file system supports hard links).**
