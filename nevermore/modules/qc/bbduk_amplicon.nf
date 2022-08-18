process qc_bbduk_stepwise_amplicon {
	label 'bbduk'

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"), emit: orphans, optional: true
    path("stats/qc/bbduk/${sample.id}/${sample.id}.*bbduk_stats.txt"), optional: true
    path("stats/qc/bbduk/${sample.id}/${sample.id}*lhist.txt"), emit: read_lengths, optional: true

    script:
	def maxmem = task.memory.toGiga()

	if (params.primers) {
		trim_params = "literal=${params.primers} minlength=${params.qc_minlen}"
	} else {
		trim_params = "ref=${adapters} minlength=${params.qc_minlen}"
	}

	def bbduk_call = "bbduk.sh -Xmx${maxmem} t=${task.cpus} ordered=t trd=t"

	ref_p5_r1 = (params.primers) ? "literal=" + params.primers.split(",")[0] : "ref=${adapters}"
	ref_p5_r2 = (params.primers && !sample.is_paired) ? "literal=" + params.primers.split(",")[1] : "ref=${adapters}"
	ref_p3_r1 = ref_p5_r2
	ref_p3_r2 = ref_p5_r1

	if (!sample.is_paired) {
		"""
	    mkdir -p ${sample.id}/
		mkdir -p stats/qc/bbduk/
		mkdir -p qc_reads/
		${bbduk_call} ${trim_params} in1=${sample.id}_R1.fastq.gz out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz stats=stats/qc/bbduk/${sample.id}/${sample.id}.fwd_bbduk_stats.txt lhist=stats/qc/bbduk/${sample.id}/${sample.id}.p5_lhist.txt
		${bbduk_call} in1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz lhist=stats/qc/bbduk/${sample.id}/${sample.id}_R1.post_lhist.txt
		"""
	} else if (params.long_reads) {
		"""
	    mkdir -p ${sample.id}/
		mkdir -p stats/qc/bbduk/
		mkdir -p qc_reads/
		mkdir -p tmp/
		${bbduk_call} ${ref_p5_r1} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R1.fastq.gz out1=fwd_p5.fastq.gz
		${bbduk_call} ${ref_p5_r2} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R2.fastq.gz out1=rev_p5.fastq.gz
		${bbduk_call} ${ref_p3_r1} minlength=${params.qc_minlen} ${params.p3_primer_params} in1=fwd_p5.fastq.gz out1=fwd.fastq.gz
		${bbduk_call} ${ref_p3_r2} minlength=${params.qc_minlen} ${params.p3_primer_params} in1=rev_p5.fastq.gz out1=rev.fastq.gz
		gzip -dc fwd.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/1//' | sort -T tmp/ > fwd.txt
        gzip -dc rev.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/2//' | sort -T tmp/ > rev.txt
		join -1 1 -2 1 fwd.txt rev.txt > both.txt
		seqtk subseq fwd.fastq.gz both.txt | gzip -c - > qc_reads/${sample.id}/${sample.id}_R1.fastq.gz
		seqtk subseq rev.fastq.gz both.txt | gzip -c - > qc_reads/${sample.id}/${sample.id}_R2.fastq.gz
		${bbduk_call} in1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz lhist=stats/qc/bbduk/${sample.id}/${sample.id}_R1.post_lhist.txt
		${bbduk_call} in1=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz lhist=stats/qc/bbduk/${sample.id}/${sample.id}_R2.post_lhist.txt
		"""
	} else {
		"""
	    mkdir -p ${sample.id}/
		mkdir -p stats/qc/bbduk/
		mkdir -p qc_reads/
		mkdir -p tmp/
		${bbduk_call} ${ref_p5_r1} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R1.fastq.gz out1=fwd.fastq.gz
		${bbduk_call} ${ref_p5_r2} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R2.fastq.gz out1=rev.fastq.gz
		gzip -dc fwd.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/1//' | sort -T tmp/ > fwd.txt
        gzip -dc rev.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/2//' | sort -T tmp/ > rev.txt
		join -1 1 -2 1 fwd.txt rev.txt > both.txt
		seqtk subseq fwd.fastq.gz both.txt | gzip -c - > qc_reads/${sample.id}/${sample.id}_R1.fastq.gz
		seqtk subseq rev.fastq.gz both.txt | gzip -c - > qc_reads/${sample.id}/${sample.id}_R2.fastq.gz
		${bbduk_call} in1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz lhist=stats/qc/bbduk/${sample.id}/${sample.id}_R1.post_lhist.txt
		${bbduk_call} in1=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz lhist=stats/qc/bbduk/${sample.id}/${sample.id}_R2.post_lhist.txt
		"""
	}
}
