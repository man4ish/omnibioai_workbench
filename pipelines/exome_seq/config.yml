process Align {
    container 'dnaseq:latest'
    input:
    file reads from params.reads
    output:
    file "*.bam" into bam_files

    """
    bwa mem ${params.ref} $reads | samtools sort -o ${reads.baseName}.bam
    samtools index ${reads.baseName}.bam
    """
}

process MarkDup {
    container 'dnaseq:latest'
    input:
    file bam from bam_files
    output:
    file "*.dedup.bam" into dedup_bam

    """
    picard MarkDuplicates I=$bam O=${bam.baseName}.dedup.bam M=${bam.baseName}.metrics.txt
    samtools index ${bam.baseName}.dedup.bam
    """
}

process BaseRecalibration {
    container 'dnaseq:latest'
    input:
    file bam from dedup_bam
    output:
    file "*.recal.bam" into recal_bam

    """
    gatk BaseRecalibrator -I $bam -R ${params.ref} -O ${bam.baseName}.recal.table
    gatk ApplyBQSR -R ${params.ref} -I $bam --bqsr-recal-file ${bam.baseName}.recal.table -O ${bam.baseName}.recal.bam
    """
}

process VariantCalling {
    container 'dnaseq:latest'
    input:
    file bam from recal_bam
    output:
    file "*.vcf" into variants

    """
    gatk HaplotypeCaller -R ${params.ref} -I $bam -O ${bam.baseName}.vcf
    """
}
