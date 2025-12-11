from integrations.igv import IGVSnapshotService

def generate_variant_view(bam, vcf, chrom, pos):
    igv = IGVSnapshotService()

    snapshot = igv.snapshot_variant(
        bam=bam,
        vcf=vcf,
        genome="hg38",
        chrom=chrom,
        start=pos-50,
        end=pos+50,
        out_name=f"{chrom}_{pos}"
    )

    return snapshot
