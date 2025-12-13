version 1.0

# Task to trim adapter sequences from FASTQ files
task trim {
    input {
        File fastq1                   # Forward FASTQ file
        File fastq2                   # Reverse FASTQ file
        File adapters                 # Adapter sequence file
        String prefix                 # Prefix for output files
        Int skewer_threads            # Number of threads for skewer
        Int minimum_read_length       # Minimum read length after trimming
        String base = basename(fastq1)  # Base name derived from forward FASTQ file
    }

    parameter_meta {
        fastq1: "Forward FASTQ file"
        fastq2: "Reverse FASTQ file"
        adapters: "Adapter sequence file"
        skewer_threads: "Number of threads for skewer"
        minimum_read_length: "Minimum read length after trimming"
        prefix: "Prefix for output files"
    }

    command {
        set -exo pipefail
        skewer -t ${skewer_threads} -y ${adapters} -l ${minimum_read_length} ${fastq1} ${fastq2} -o ${prefix}
    }

    output {
        File out_pair1 = '${prefix}-trimmed-pair1.fastq'  # Trimmed forward FASTQ file
        File out_pair2 = '${prefix}-trimmed-pair2.fastq'  # Trimmed reverse FASTQ file
    }

    runtime {
        # Docker image for task execution
        docker: "docker.io/man4ish/rnaseq:latest"
    }
}

# Task for quantification using Kallisto
task quantification {
    input {
        File fastq1                   # Trimmed forward FASTQ file
        File fastq2                   # Trimmed reverse FASTQ file
        File idx                      # Index file for Kallisto
        File gtf                      # GTF annotation file
        Int kallisto_threads          # Number of threads for Kallisto
        Int bootstrap_samples         # Number of bootstrap samples
        String base = basename(fastq1)  # Base name derived from forward FASTQ file
    }

    parameter_meta {
        fastq1: "Trimmed forward FASTQ file"
        fastq2: "Trimmed reverse FASTQ file"
        idx: "Index file for Kallisto"
        gtf: "GTF annotation file"
        kallisto_threads: "Number of threads for Kallisto"
        bootstrap_samples: "Number of bootstrap samples"
    }

    command {
        set -exo pipefail
        mkdir ${base}_output
        /software/utils/kallisto quant -t ${kallisto_threads} -b ${bootstrap_samples} --rf-stranded --genomebam -i ${idx} -g ${gtf} ${fastq1} ${fastq2} -o ${base}_output
        tar -cvzf ${base}_kallisto_output.tar.gz ${base}_output
    }

    output {
        File out = '${base}_kallisto_output.tar.gz'  # Kallisto output archive
    }

    runtime {
        # Docker image for task execution
        docker: "docker.io/man4ish/rnaseq:latest"
    }
}

# Task to align reads using STAR
task align {
    input {
        File ref_tar                   # Reference genome archive
        File fastq1                   # Trimmed forward FASTQ file
        File fastq2                   # Trimmed reverse FASTQ file
        Int STAR_threads              # Number of threads for STAR
        String base = basename(ref_tar, ".tar.gz")  # Base name derived from reference genome file
    }

    parameter_meta {
        fastq1: "Trimmed forward FASTQ file"
        fastq2: "Trimmed reverse FASTQ file"
        ref_tar: "Reference genome archive"
        STAR_threads: "Number of threads for STAR"
    }

    command {
        set -exo pipefail
        mkdir ref_bundle
        tar -xzf ${ref_tar} -C ref_bundle --no-same-owner
        ref_path=$(realpath ref_bundle)
        mv ref_bundle/*/* ref_bundle
        /software/utils/STAR --genomeDir ${ref_path} --runThreadN ${STAR_threads} --outSAMtype BAM Unsorted --readFilesIn ${fastq1} ${fastq2} --outFileNamePrefix ${base}_sample
        cp ${base}_sampleAligned.out.bam ${base}_sample.bam
    }

    output {
        File out_bam = '${base}_sample.bam'  # Aligned BAM file
    }

    runtime {
        # Docker image for task execution
        docker: "docker.io/man4ish/rnaseq:latest"
    }
}

# Task to sort and index BAM files
task sort_index {
    input {
        File bam_file                  # BAM file to be sorted and indexed
        String base = basename(bam_file)  # Base name derived from BAM file
    }

    parameter_meta {
        bam_file: "BAM file to be sorted and indexed"
    }

    command {
        set -exo pipefail
        samtools sort ${bam_file} > sorted_${base}
        samtools index sorted_${base}
    }

    output {
        File out_sortedbam = 'sorted_${base}'  # Sorted BAM file
        File out_sortedbam_index = 'sorted_${base}.bai'  # Index file for sorted BAM
    }

    runtime {
        # Docker image for task execution
        docker: "docker.io/man4ish/rnaseq:latest"
    }
}

# Task to generate summary metrics using Picard
task gen_summary {
    input {
        File bam_file                  # Sorted BAM file
        File ref_flat                  # Reference flat file for Picard
        File ribosomal_interval        # Ribosomal interval file for Picard
        File ref_seq                   # Reference sequence file
        String base = basename(bam_file)  # Base name derived from BAM file
    }

    parameter_meta {
        bam_file: "Sorted BAM file"
        ref_flat: "Reference flat file for Picard"
        ribosomal_interval: "Ribosomal interval file for Picard"
        ref_seq: "Reference sequence file"
    }

    command {
        set -exo pipefail
        java -jar /software/utils/picard.jar CollectRnaSeqMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS \
            REF_FLAT=${ref_flat} RIBOSOMAL_INTERVALS=${ribosomal_interval} \
            STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
            CHART_OUTPUT=${base}_position.vs.coverage.plot.pdf \
            INPUT=${bam_file} OUTPUT=${base}.rna.summary \
            REFERENCE_SEQUENCE=${ref_seq} VALIDATION_STRINGENCY=LENIENT
    }

    output {
        File out_summary = '${base}.rna.summary'  # RNA-seq summary file
        File out_plots = '${base}_position.vs.coverage.plot.pdf'  # Coverage plot
    }

    runtime {
        # Docker image for task execution
        docker: "docker.io/man4ish/rnaseq:latest"
    }
}

# Main workflow for RNA-seq analysis
workflow Rnaseq {
    input {
        File fastq1                   # Forward FASTQ file
        File fastq2                   # Reverse FASTQ file
        File adapters                 # Adapter sequence file
        Int skewer_threads            # Number of threads for skewer
        Int minimum_read_length       # Minimum read length after trimming
        String prefix                 # Prefix for output files

        Int kallisto_threads          # Number of threads for Kallisto
        Int bootstrap_samples         # Number of bootstrap samples
        File idx                      # Index file for Kallisto
        File gtf                      # GTF annotation file

        Int STAR_threads              # Number of threads for STAR
        File ref_tar                  # Reference genome archive

        File ref_flat                 # Reference flat file for Picard
        File ribosomal_interval       # Ribosomal interval file for Picard
        File ref_seq                  # Reference sequence file
    }

    call trim {
        input: skewer_threads=skewer_threads, minimum_read_length=minimum_read_length, fastq1=fastq1, fastq2=fastq2, adapters=adapters, prefix=prefix
    }

    call quantification {
        input: fastq1=trim.out_pair1, fastq2=trim.out_pair2, kallisto_threads=kallisto_threads, bootstrap_samples=bootstrap_samples, idx=idx, gtf=gtf
    }

    call align {
        input: fastq1=trim.out_pair1, fastq2=trim.out_pair2, STAR_threads=STAR_threads, ref_tar=ref_tar
    }

    call sort_index {
        input: bam_file=align.out_bam
    }

    call gen_summary {
        input: bam_file=sort_index.out_sortedbam, ref_flat=ref_flat, ribosomal_interval=ribosomal_interval, ref_seq=ref_seq
    }
}