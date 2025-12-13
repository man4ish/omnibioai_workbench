#!/usr/bin/env nextflow

params.fastq_dir = "data/sc_rnaseq/"
params.genome = "refdata-cellranger-GRCh38-3.0.0"
params.outdir = "results/sc_rnaseq"

process CellRangerCount {
    container 'scrnaseq:latest'
    input:
    tuple val(sample_id), file(fastqs)

    output:
    file("${sample_id}") into sc_matrices

    """
    cellranger count --id=${sample_id} \
                     --transcriptome=${params.genome} \
                     --fastqs=${fastqs} \
                     --sample=${sample_id} \
                     --localcores 16 --localmem 64
    """
}

process ScanpyQC {
    container 'scrnaseq:latest'
    input:
    file matrix from sc_matrices
    output:
    file "*.h5ad" into h5ad_files

    """
    python - <<EOF
import scanpy as sc
adata = sc.read_10x_h5("${matrix}/outs/filtered_feature_bc_matrix.h5")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.write("${matrix}.h5ad")
EOF
    """
}

workflow {
    Channel.fromPath("${params.fastq_dir}/*").map { f -> tuple(f.baseName, f) } | CellRangerCount()
    ScanpyQC()
}
