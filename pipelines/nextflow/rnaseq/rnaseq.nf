#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq1
params.fastq2
params.adapters
params.prefix
params.skewer_threads = 4
params.minimum_read_length = 30

params.idx
params.gtf
params.kallisto_threads = 4
params.bootstrap_samples = 100

params.ref_tar
params.STAR_threads = 8

params.ref_flat
params.ribosomal_interval
params.ref_seq

process TRIM {
    tag "$params.prefix"

    container 'docker.io/man4ish/rnaseq:latest'

    input:
    path fastq1 from params.fastq1
    path fastq2 from params.fastq2
    path adapters from params.adapters
    val prefix from params.prefix
    val threads from params.skewer_threads
    val min_len from params.minimum_read_length

    output:
    path "${prefix}-trimmed-pair1.fastq", emit: trimmed1
    path "${prefix}-trimmed-pair2.fastq", emit: trimmed2

    script:
    """
    skewer -t ${threads} -y ${adapters} -l ${min_len} ${fastq1} ${fastq2} -o ${prefix}
    """
}

process KALLISTO_QUANT {
    container 'docker.io/man4ish/rnaseq:latest'

    input:
    path fastq1 from TRIM.trimmed1
    path fastq2 from TRIM.trimmed2
    path idx from params.idx
    path gtf from params.gtf
    val threads from params.kallisto_threads
    val bootstrap from params.bootstrap_samples
    val prefix from params.prefix

    output:
    path "${prefix}_kallisto_output.tar.gz", emit: kallisto_out

    script:
    """
    mkdir ${prefix}_output
    /software/utils/kallisto quant -t ${threads} -b ${bootstrap} --rf-stranded --genomebam -i ${idx} -g ${gtf} ${fastq1} ${fastq2} -o ${prefix}_output
    tar -cvzf ${prefix}_kallisto_output.tar.gz ${prefix}_output
    """
}

process STAR_ALIGN {
    container 'docker.io/man4ish/rnaseq:latest'

    input:
    path fastq1 from TRIM.trimmed1
    path fastq2 from TRIM.trimmed2
    path ref_tar from params.ref_tar
    val threads from params.STAR_threads
    val prefix from params.prefix

    output:
    path "${prefix}_sample.bam", emit: aligned_bam

    script:
    """
    mkdir ref_bundle
    tar -xzf ${ref_tar} -C ref_bundle --no-same-owner
    ref_path=\$(realpath ref_bundle)
    mv ref_bundle/*/* ref_bundle
    /software/utils/STAR --genomeDir \${ref_path} --runThreadN ${threads} --outSAMtype BAM Unsorted --readFilesIn ${fastq1} ${fastq2} --outFileNamePrefix ${prefix}_sample
    cp ${prefix}_sampleAligned.out.bam ${prefix}_sample.bam
    """
}

process SORT_INDEX {
    container 'docker.io/man4ish/rnaseq:latest'

    input:
    path bam_file from STAR_ALIGN.aligned_bam
    val prefix from params.prefix

    output:
    path "sorted_${prefix}_sample.bam", emit: sorted_bam
    path "sorted_${prefix}_sample.bam.bai", emit: bam_index

    script:
    """
    samtools sort ${bam_file} > sorted_${prefix}_sample.bam
    samtools index sorted_${prefix}_sample.bam
    """
}

process PICARD_SUMMARY {
    container 'docker.io/man4ish/rnaseq:latest'

    input:
    path bam_file from SORT_INDEX.sorted_bam
    path ref_flat from params.ref_flat
    path ribosomal_interval from params.ribosomal_interval
    path ref_seq from params.ref_seq
    val prefix from params.prefix

    output:
    path "${prefix}.rna.summary"
    path "${prefix}_position.vs.coverage.plot.pdf"

    script:
    """
    java -jar /software/utils/picard.jar CollectRnaSeqMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS \
        REF_FLAT=${ref_flat} RIBOSOMAL_INTERVALS=${ribosomal_interval} \
        STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
        CHART_OUTPUT=${prefix}_position.vs.coverage.plot.pdf \
        INPUT=${bam_file} OUTPUT=${prefix}.rna.summary \
        REFERENCE_SEQUENCE=${ref_seq} VALIDATION_STRINGENCY=LENIENT
    """
}

workflow {
    TRIM()
    KALLISTO_QUANT()
    STAR_ALIGN()
    SORT_INDEX()
    PICARD_SUMMARY()
}
