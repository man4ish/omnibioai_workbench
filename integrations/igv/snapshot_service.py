from pathlib import Path
from .client import IGVClient

class IGVSnapshotService:
    """
    Automates:
      - Load tracks
      - Jump to loci
      - Take snapshots
    """

    def __init__(self, snapshot_dir="data/igv_snapshots"):
        self.client = IGVClient()
        self.snapshot_dir = Path(snapshot_dir)
        self.snapshot_dir.mkdir(parents=True, exist_ok=True)

        # Configure IGV snapshot directory
        self.client.snapshot_directory(str(self.snapshot_dir))

    def snapshot_variant(self, bam, vcf, genome, chrom, start, end, out_name):
        self.client.set_genome(genome)
        self.client.load(bam)
        self.client.load(vcf)

        region = f"{chrom}:{start}-{end}"
        self.client.goto(region)

        fname = f"{out_name}.png"
        self.client.snapshot(fname)

        return str(self.snapshot_dir / fname)
