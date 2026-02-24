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

