# Usage

## CLI
```
nextflow run /path/to/vortex_knight --input_dir </path/to/read_files> --output_dir </path/to/output_dir> [PARAMETERS]
```

## Input data

* Fastq files need to be ordered into a sample-per-folder structure (s. tree below).
* Supported file endings are `.fastq,.fq,.fastq.gz,.fq.gz,.fastq.bz2,.fq.bz2`
* Files in sample-specific folders will automatically be assigned to a sample labeled with the folder name.
* Filenames for paired-end data need to share a common prefix terminated by `_R?[12]` in order to be automatically matched. 
* Suffixes following the `_R?[12]` pattern, such as `_001` from newer Illumina machines need to be removed.
* In case of paired-end data with additional unpaired reads (e.g. as obtained from preprocessed ENA datasets), the unpaired files can be automatically picked up from the same folder but should be labeled as `<sample_prefix>.single(s).<fastq-suffix>`. 

In the following example, the path to the input dataset needs to be set as `--input_dir /path/to/read_files`. `vortex_knight` will then automatically detect the sample-specific folders within that directory.

```
/path/to/read_files
└── sample1
  ├── sample1_R1.fastq.gz
  └── sample2_R2.fastq.gz
```

## Parameters
* `run_preprocessing [true]`: Run preprocessing (quality control, [human] host removal, rRNA removal).
* `remove_host [true]`: Run host removal.
* `drop_orphans [false]`: Drop paired-end reads whose mate did not survive quality processing.
* `qc_minlen [45]`: Drop reads shorter than `qc_minlen` base pairs.
* `qc_params_shotgun ["qtrim=rl trimq=3 maq=25 ktrim=r k=23 mink=11 hdist=1 ftm=5 entropy=0.5 entropywindow=50 entropyk=5 tpe tbo"]`: bbduk parameter string.
* `remove_host_kraken2_db`: Path to a kraken2 database for host removal.
* `kraken2_min_hit_groups [10]`: kraken2 sensitivity cutoff.


### Input parameters

* `--input_dir` should be a folder with bam files or with gzipped fastq files. For fastq files, individual samples should be separated into individual folders.
* `--output_dir` is `vknight_out` in the local directory by default.
* `--skip_<analysis>`, `--run_<analysis>` skips, resp. explicitly requires execution of the specified analysis (`motus`, `pathseq`, `count_reads`, `mtags`, `mapseq`, `kraken2`)
* `--publishMode` allows to switch between various modes of how results files are placed in the `output_dir` (cf. NextFlow documentation)

#### Notes
* `mapseq` can only run in combination with `mtags` and when the parameter `mapseq_bin` is explicitly set.
* `kraken2` can only run when the parameter `kraken_database` is set.
* `pathseq` can only run when the parameter `pathseq_database` is set.
* a pre-downloaded motus database can be set with the parameter `motus_database`.
* results are only collated if the parameter `collate_script` is set. (TODO -> change to baseDir?)


