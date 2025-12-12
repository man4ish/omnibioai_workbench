# tests/integrations/test_igv_plot.py
import os
from pathlib import Path
from integrations.igv.client import IGVClient
from integrations.igv.exceptions import IGVConnectionError, IGVCommandError

def plot_igv_example():
    client = IGVClient(host="localhost", port=60151, timeout=5)

    # Directory to save snapshots
    snapshot_dir = Path("data/igv_snapshots")
    snapshot_dir.mkdir(exist_ok=True)

    try:
        # 1. Set genome
        client.set_genome("hg38")

        # 2. Load a BED file (genomic regions)
        bed_file = "tests/data/example.bed"  # replace with your BED/BAM/VCF path
        client.load(bed_file)

        # 3. Go to a specific region
        client.goto("chr1:1000000-1010000")

        # 4. Set snapshot directory
        client.snapshot_directory(str(snapshot_dir))

        # 5. Take a snapshot
        snapshot_file = snapshot_dir / "igv_region.png"
        client.snapshot(str(snapshot_file))
        print(f"Snapshot saved at {snapshot_file}")

    except IGVConnectionError:
        print("IGV is not running or reachable.")
    except IGVCommandError as e:
        print(f"IGV command failed: {e}")

if __name__ == "__main__":
    plot_igv_example()

