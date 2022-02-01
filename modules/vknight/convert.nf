process bam2fq {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("fastq/${sample}/${sample}*.fastq.gz"), emit: reads

    script:
    """
    set -o pipefail
    mkdir -p fastq/${sample}
    samtools collate -@ $task.cpus -u -O $bam | samtools fastq -F 0x900 -0 ${sample}_other.fastq.gz -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz

    if [[ "\$?" -eq 0 ]];
    then

        if [[ -z "\$(gzip -dc ${sample}_R1.fastq.gz | head -n 1)" ]];
        then
            if [[ ! -z "\$(gzip -dc ${sample}_other.fastq.gz | head -n 1)" ]];
            then
                mv -v ${sample}_other.fastq.gz fastq/${sample}/${sample}_R1.fastq.gz;
            fi;
        else
                mv -v ${sample}_R1.fastq.gz fastq/${sample}/;
                if [[ ! -z "\$(gzip -dc ${sample}_R2.fastq.gz | head -n 1)" ]];
                then
                    mv -v ${sample}_R2.fastq.gz fastq/${sample}/;
                fi;
        fi;

        ls -l *.fastq.gz
        ls -l fastq/${sample}/*.fastq.gz
        rm -rf *.fastq.gz
    fi;
    """
}

process fq2bam {
    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("out/${sample.id}.bam"), emit: reads

    script:
    if (sample.is_paired) {
        """
        mkdir -p out
        gatk FastqToSam -F1 ${fq[0]} -F2 ${fq[1]} -O out/${sample.id}.bam -SM ${sample.id}
        """
    } else {
        """
        mkdir -p out
        gatk FastqToSam -F1 ${fq[0]} -O out/${sample.id}.bam -SM ${sample.id}
        """
    }
}

process prepare_fastqs {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("fastq/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads

    script:
    if (sample.is_paired) {
        """
        mkdir -p fastq/${sample.id}
        ln -sf ../../${fq[0]} fastq/${sample.id}/${sample.id}_R1.fastq.gz
        ln -sf ../../${fq[1]} fastq/${sample.id}/${sample.id}_R2.fastq.gz
        """
    } else {
        """
        mkdir -p fastq/${sample.id}
        ln -sf ../../${fq[0]} fastq/${sample.id}/${sample.id}_R1.fastq.gz
        """
    }
}
