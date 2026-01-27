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
// - copies the leading "#mpa_..." line from the first input counts file (if present)
// - writes header as: clade_name <sample>.mp4 ...
// - reads per-sample counts files that look like:
//     #mpa_...
//     clade_name    estimated_number_of_reads_from_the_clade
//     k__Bacteria   48678090
////////////////////////////////////////////////////////////////////////////////
process collate_metaphlan4_counts {
  publishDir params.output_dir, mode: "copy"
  container "quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
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

  files = "${count_tables}".split()
  if not files:
    raise SystemExit("ERROR: No count tables provided")

  # Copy the MetaPhlAn database/version header if present
  mpa_line = None
  with open(files[0], "r") as f0:
    for line in f0:
      if line.startswith("#mpa_"):
        mpa_line = line.rstrip("\\n")
        break

  def sample_from_path(p):
    base = os.path.basename(p)
    s = re.sub(r"\\.mp4\\.counts\\.tsv\$", "", base)
    return f"{s}.mp4"

  counts = defaultdict(dict)  # clade -> {sample: value}
  samples = []

  for fp in files:
    s = sample_from_path(fp)
    samples.append(s)

    with open(fp, "r", newline="") as f:
      # Skip initial comment lines (e.g. #mpa_...)
      first_data_line = None
      for line in f:
        if line.startswith("#"):
          continue
        first_data_line = line
        break

      if first_data_line is None:
        raise SystemExit(f"ERROR: No header/data found in {fp}")

      # Build a DictReader starting from the header line we just found
      header = first_data_line.rstrip("\\n").split("\\t")

      # Determine column names (support both old/new extractors)
      # Preferred:
      clade_col = "clade_name" if "clade_name" in header else ("clade" if "clade" in header else None)
      val_col = (
        "estimated_number_of_reads_from_the_clade"
        if "estimated_number_of_reads_from_the_clade" in header
        else ("estimated_reads" if "estimated_reads" in header else None)
      )

      if clade_col is None or val_col is None:
        raise SystemExit(f"ERROR: Unexpected header in {fp}: {header}")

      # Now read the rest of the file as TSV with that header
      reader = csv.DictReader(f, fieldnames=header, delimiter="\\t")
      for row in reader:
        if not row:
          continue
        clade = row.get(clade_col)
        val = row.get(val_col)
        if clade is None or val is None:
          continue
        counts[clade][s] = val

  samples = sorted(set(samples))

  with open("metaphlan4_counts_matrix.tsv", "w", newline="") as out:
    if mpa_line:
      out.write(mpa_line + "\\n")
    w = csv.writer(out, delimiter="\\t")
    w.writerow(["clade_name"] + samples)
    for clade in sorted(counts.keys()):
      w.writerow([clade] + [counts[clade].get(s, "0") for s in samples])
  PY
  """
}
