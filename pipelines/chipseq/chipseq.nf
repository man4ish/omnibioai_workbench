process PeakCalling {
    container 'chipseq:latest'
    input:
    file bam from aligned_bam
    output:
    file "*.narrowPeak" into peaks

    """
    macs2 callpeak -t $bam -f BAM -g hs -n ${bam.baseName} --outdir .
    """
}

process SignalTracks {
    container 'chipseq:latest'
    input:
    file bam from aligned_bam
    output:
    file "*.bw" into bigwigs

    """
    bamCoverage -b $bam -o ${bam.baseName}.bw
    """
}
