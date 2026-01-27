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
// - preserves the leading "#mpa_..." line (if present)
// - outputs a small 2-column TSV:
//     #mpa_...
//     clade_name  estimated_number_of_reads_from_the_clade
////////////////////////////////////////////////////////////////////////////////
process extract_mp4_counts {
  publishDir params.output_dir, mode: "copy"
  container "quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
  label "metaphlan4"
  label "mini"

  input:
  tuple val(sample), path(mp4_table)

  output:
  tuple val(sample), path("${sample.id}.mp4.counts.tsv"), emit: mp4_counts

  script:
  def in_table = mp4_table
  def out_counts = "${sample.id}.mp4.counts.tsv"

  """
  python - <<'PY'
  import re

  infile = "${in_table}"
  outfile = "${out_counts}"

  mpa_line = None
  header = None
  clade_idx = None
  target_idx = None
  out = None
  wrote_any = False

  with open(infile, "r") as f:
    for line in f:
      line = line.rstrip("\\n")
      if not line:
        continue

      # Capture DB/version line (optional)
      if line.startswith("#mpa_") and mpa_line is None:
        mpa_line = line
        continue

      # Capture header row (MetaPhlAn prints it as a comment)
      if line.startswith("#clade_name"):
        header = line.lstrip("#").split("\\t")

        # Determine indices once
        clade_idx = header.index("clade_name") if "clade_name" in header else 0

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

        continue

      # Skip any other comment lines
      if line.startswith("#"):
        continue

      # First data row requires that we already parsed the header
      if header is None or clade_idx is None or target_idx is None:
        raise SystemExit(f"ERROR: Could not find '#clade_name' header in {infile}")

      # Open output on first data line (keeps behaviour stable)
      if out is None:
        out = open(outfile, "w")
        if mpa_line:
          out.write(mpa_line + "\\n")
        out.write("clade_name\\testimated_number_of_reads_from_the_clade\\n")

      parts = line.split("\\t")
      if len(parts) <= max(clade_idx, target_idx):
        continue

      out.write(parts[clade_idx] + "\\t" + parts[target_idx] + "\\n")
      wrote_any = True

  if out is not None:
    out.close()

  if not wrote_any:
    raise SystemExit(f"ERROR: No data rows found in {infile}")
  PY
  """
}


////////////////////////////////////////////////////////////////////////////////
// NEW: Merge all per-sample counts tables into one matrix (clade x sample)
// - copies the leading "#mpa_..." line from the first input counts file (if present)
// - writes header as: clade_name <sample>.mp4 ...
// - supports both column-header variants:
//     clade_name + estimated_number_of_reads_from_the_clade   (preferred)
//     clade      + estimated_reads                            (legacy)
////////////////////////////////////////////////////////////////////////////////
process collate_metaphlan4_counts {
  publishDir params.output_dir, mode: "copy"
  container "quay.io/biocontainers/metaphlan:4.1.0--pyhca03a8a_0"
  label "metaphlan4"
  label "mini"

  input:
  path(count_tables)

  output:
  path("metaphlan4_counts_matrix.tsv"), emit: mp4_counts_matrix

  script:
  def out_matrix = "metaphlan4_counts_matrix.tsv"

  """
  python - <<'PY'
  import os, re, csv
  from collections import defaultdict

  files = "${count_tables}".split()
  if not files:
    raise SystemExit("ERROR: No count tables provided")

  # Copy MetaPhlAn database/version line from the first file (optional)
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

  counts = defaultdict(dict)
  samples = []

  for fp in files:
    s = sample_from_path(fp)
    samples.append(s)

    with open(fp, "r", newline="") as f:
      # Find the first non-comment line (= header)
      header_line = None
      for line in f:
        if line.startswith("#") or not line.strip():
          continue
        header_line = line
        break

      if header_line is None:
        raise SystemExit(f"ERROR: No header/data found in {fp}")

      header = header_line.rstrip("\\n").split("\\t")

      clade_col = "clade_name" if "clade_name" in header else ("clade" if "clade" in header else None)
      val_col = (
        "estimated_number_of_reads_from_the_clade"
        if "estimated_number_of_reads_from_the_clade" in header
        else ("estimated_reads" if "estimated_reads" in header else None)
      )

      if clade_col is None or val_col is None:
        raise SystemExit(f"ERROR: Unexpected header in {fp}: {header}")

      reader = csv.DictReader(f, fieldnames=header, delimiter="\\t")
      for row in reader:
        clade = row.get(clade_col)
        val = row.get(val_col)
        if not clade or val is None:
          continue
        counts[clade][s] = val

  samples = sorted(set(samples))

  with open("${out_matrix}", "w", newline="") as out:
    if mpa_line:
      out.write(mpa_line + "\\n")
    w = csv.writer(out, delimiter="\\t")
    w.writerow(["clade_name"] + samples)
    for clade in sorted(counts.keys()):
      w.writerow([clade] + [counts[clade].get(s, "0") for s in samples])
  PY
  """
}
