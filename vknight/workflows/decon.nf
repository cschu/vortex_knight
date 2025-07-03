include { starmap } from "../modules/decon/starmap"
include { bwa_mem_align } from "../../nevermore/modules/align/bwa"
include { hostile } from "../../nevermore/modules/decon/hostile"
include { bam2fq } from "../../nevermore/modules/converters/bam2fq"
include { flagstats } from "../../nevermore/modules/stats"
include { fastqc } from "../../nevermore/modules/qc/fastqc"

// params.remove_host_star_db = "/scratch/fspringe/Databases/STAR_index/Ahmad_crc-eco"
// params.remove_host_bwa_index = "/g/scb/zeller/fspringe/ReferenceGenomes/T2T-CHM13v2.0/chm13v2.0.fa.idx"

params.decon_mode = "bwa_mem"

workflow vk_decon {

	take:
		reads_ch
	main:
		flagstats_reads_ch = Channel.empty()

		if (params.library_type == "rna") {

			starmap(reads_ch, params.remove_host_star_db)

			bam2fq(starmap.out.bam, true)

			reads_ch = bam2fq.out.reads

			flagstats_reads_ch = flagstats_reads_ch.mix(starmap.out.bam)

		} else {

			if (params.decon_mode == "bwa_mem") {

				bwa_mem_align(reads_ch, params.remove_host_bwa_index, false)

				bam2fq(bwa_mem_align.out.bam, true)

				reads_ch = bam2fq.out.reads

				flagstats_reads_ch = flagstats_reads_ch.mix(bwa_mem_align.out.bam)

			} else if (params.decon_mode == "hostile") {

				hostile(reads_ch, params.hostile_db)

				reads_ch = hostile.out.reads

				// flagstats?

			}


		}

		flagstats(flagstats_reads_ch, "decon")
		fastqc(reads_ch, "decon_survivors")

	emit:
		reads = reads_ch
		


}