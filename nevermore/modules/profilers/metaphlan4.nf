process run_metaphlan4 {
  publishDir params.output_dir, mode: "copy"
  container "quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
  tag "${sample.id}"
  label "process_high"
  label "metaphlan4"

  input:
  tuple val(sample), path(fastqs)
  path(mp4_db)

  output:
  tuple val(sample), path("${sample.id}.mp4.txt"), emit: mp4_table
  tuple val(sample), path("${sample.id}.mp4.sam.bz2"), emit: mp4_sam, optional: (params.run_samestr || params.samestr_compatible_output)
  // tuple val(sample), path("${sample.id}.bowtie2.bz2"), emit: mp4_bt2

  script:
  def mp4_params = "--bowtie2db ${mp4_db} --input_type fastq --nproc ${task.cpus} --tmp_dir tmp/"
  def mp4_input = ""
  def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"

  def samestr_params = ""
  if (params.run_samestr || params.samestr_compatible_output) {
    samestr_params = "--samout ${sample.id}.mp4.sam.bz2"
  }

  def input_files = []
  input_files += fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
  input_files += fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
  input_files += fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
  mp4_input = input_files.join(',')

  """
  mkdir -p tmp/

  # changed: rel_ab -> rel_ab_w_read_stats (adds estimated reads per clade)
  metaphlan ${mp4_input} ${mp4_params} ${bt2_out} -o ${sample.id}.mp4.txt ${samestr_params} -t rel_ab_w_read_stats

  # keep your behaviour (creates file even when no samout was requested)
  touch ${sample.id}.mp4.sam.bz2
  """
}


process combine_metaphlan4 {
  container "quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
  label "metaphlan4"

  input:
  tuple val(sample), path(bt2)

  output:
  tuple val(sample), path("metaphlan4/${sample.id}/${sample.id}.mp4.txt"), emit: mp4_table

  script:
  def mp4_params = "--input_type bowtie2out --nproc ${task.cpus} --tmp_dir tmp/"
  def mp4_input = "${sample.id}.bowtie2.bz2,${sample.id}.singles.bowtie2.bz2"
  """
  mkdir -p metaphlan4/${sample.id}/

  metaphlan ${mp4_input} ${mp4_params} -o metaphlan4/${sample.id}/${sample.id}.mp4.txt
  """
}


process collate_metaphlan4_tables {
  publishDir params.output_dir, mode: "copy"
  container "quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
  label "metaphlan4"
  label "mini"

  input:
  path(tables)

  output:
  path("metaphlan4_abundance_table.txt"), emit: mp4_abundance_table

  script:
  """
  mkdir -p metaphlan4/

  merge_metaphlan_tables.py ${tables} > metaphlan4_abundance_table.txt
  """
}


////////////////////////////////////////////////////////////////////////////////
// NEW: Extract per-sample estimated-read "counts" from rel_ab_w_read_stats table
////////////////////////////////////////////////////////////////////////////////
process extract_mp4_counts {
  publishDir params.output_dir, mode: params.publish_mode ?: "copy"
  container "quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
  label "metaphlan4"
  label "mini"

  input:
  tuple val(sample), path(mp4_table)

  output:
  tuple val(sample), path("${sample.id}.mp4.counts.tsv"), emit: mp4_counts

  script:
  """
  python - <<'PY'
  import re

  infile = "${mp4_table}"
  sample = "${sample.id}"
  outfile = f"{sample}.mp4.counts.tsv"

  header = None
  data_started = False

  with open(infile, "r") as f:
    for line in f:
      line = line.rstrip("\\n")
      if not line:
        continue

      # MetaPhlAn header is a comment line starting with '#clade_name'
      if line.startswith("#clade_name"):
        header = line.lstrip("#").split("\\t")
        continue

      # Skip other comments
      if line.startswith("#"):
        continue

      # First non-comment data line begins after we have header
      if header is None:
        # If file is in unexpected format, fail early with a helpful message
        raise SystemExit(f"ERROR: Could not find '#clade_name' header in {infile}")

      # Now we're reading data rows
      if not data_started:
        # Determine column indices once
        clade_idx = header.index("clade_name") if "clade_name" in header else 0

        # Find estimated reads column by exact name (MetaPhlAn 4)
        # fallback to regex if name changes slightly
        if "estimated_number_of_reads_from_the_clade" in header:
          target_idx = header.index("estimated_number_of_reads_from_the_clade")
        else:
          target_idx = None
          for j, col in enumerate(header):
            if re.search(r"estimated.*reads.*clade", col, re.IGNORECASE):
              target_idx = j
              break
          if target_idx is None:
            raise SystemExit("ERROR: Could not find estimated reads column in header: " + str(header))

        out = open(outfile, "w")
        out.write("clade\\testimated_reads\\n")
        data_started = True

      parts = line.split("\\t")
      if len(parts) <= max(clade_idx, target_idx):
        continue
      out.write(parts[clade_idx] + "\\t" + parts[target_idx] + "\\n")

  if not data_started:
    raise SystemExit(f"ERROR: No data rows found in {infile}")

  out.close()
  PY
  """
}



////////////////////////////////////////////////////////////////////////////////
// NEW: Merge all per-sample counts tables into one matrix (clade x sample)
////////////////////////////////////////////////////////////////////////////////
process collate_metaphlan4_counts {
  publishDir params.output_dir, mode: "copy"
  container "python:3.11-slim"
  label "mini"

  input:
  path(count_tables)

  output:
  path("metaphlan4_counts_matrix.tsv"), emit: mp4_counts_matrix

  script:
  """
  python - <<'PY'
  import os, re, csv
  from collections import defaultdict

  # Nextflow passes a whitespace-separated list of paths in the interpolation below
  files = "${count_tables}".split()
  if not files:
    raise SystemExit("ERROR: No count tables provided")

  def sample_from_path(p):
    base = os.path.basename(p)
    return re.sub(r"\\.mp4\\.counts\\.tsv\$", "", base)

  counts = defaultdict(dict)  # clade -> {sample: value}
  samples = []

  for fp in files:
    s = sample_from_path(fp)
    samples.append(s)
    with open(fp, newline="") as f:
      rdr = csv.DictReader(f, delimiter="\\t")
      for row in rdr:
        clade = row["clade"]
        val = row["estimated_reads"]
        counts[clade][s] = val

  samples = sorted(set(samples))

  with open("metaphlan4_counts_matrix.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\\t")
    w.writerow(["clade"] + samples)
    for clade in sorted(counts.keys()):
      w.writerow([clade] + [counts[clade].get(s, "0") for s in samples])
  PY
  """
}

