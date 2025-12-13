#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "data/rnaseq/*_R{1,2}.fastq.gz"
params.genome_index = "genomes/hg38_star_index"
params.gtf = "genomes/hg38.gtf"
params.outdir = "results/rnaseq"

process FastQC {
    container 'rnaseq:latest'
    input:
    file reads from params.reads
    output:
    file "*.html" into qc_reports

    """
    fastqc $reads -o .
    """
}

process TrimAdapters {
    container 'rnaseq:latest'
    input:
    file reads from params.reads
    output:
    file "*.trimmed.fastq.gz" into trimmed_reads

    """
    trim_galore --paired $reads
    """
}

process AlignSTAR {
    container 'rnaseq:latest'
    input:
    file reads from trimmed_reads
    output:
    file "*.bam" into aligned_bam

    """
    STAR --genomeDir ${params.genome_index} \
         --readFilesIn $reads \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${reads.baseName}.
    """
}

process QuantifyGenes {
    container 'rnaseq:latest'
    input:
    file bam from aligned_bam
    output:
    file "counts.txt" into gene_counts

    """
    featureCounts -a ${params.gtf} -o counts.txt $bam
    """
}

process DEAnalysis {
    container 'rnaseq:latest'
    input:
    file counts from gene_counts.collect()
    output:
    file "DE_results.csv"

    """
    Rscript -e '
      library(DESeq2);
      counts <- read.table("counts.txt", header=T, row.names=1);
      coldata <- data.frame(row.names=colnames(counts), condition=c("A","B"));
      dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition);
      dds <- DESeq(dds);
      res <- results(dds);
      write.csv(as.data.frame(res), file="DE_results.csv")
    '
    """
}

workflow {
    FastQC()
    TrimAdapters()
    AlignSTAR()
    QuantifyGenes()
    DEAnalysis()
}
