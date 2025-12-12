import os
from pathlib import Path
from integrations.igv.client import IGVClient
from integrations.igv.exceptions import IGVConnectionError, IGVCommandError
from omnibioai.services.network_viz import NetworkViz

# -----------------------------
# 1. Setup IGV client
# -----------------------------
igv = IGVClient(host="localhost", port=60151)

try:
    sock = igv._connect()
    sock.close()
except IGVConnectionError:
    print("IGV is not running or socket server not enabled. Please start IGV with Tools → Advanced → Enable Socket Server.")
    exit(1)

# -----------------------------
# 2. Set genome in IGV
# -----------------------------
genome = "hg38"
resp = igv.set_genome(genome)
print(f"Set genome response: {resp}")

# -----------------------------
# 3. Create a temporary BED track
# -----------------------------
tmp_dir = Path("tmp_igv_test")
tmp_dir.mkdir(exist_ok=True)
bed_file = tmp_dir / "test_track.bed"
bed_file.write_text(
    "chr1\t10000\t10100\tfeature1\n"
    "chr1\t20000\t20100\tfeature2\n"
    "chr2\t30000\t30100\tfeature3\n"
)

resp = igv.load(str(bed_file))
print(f"Load track response: {resp}")

# -----------------------------
# 4. Go to locus
# -----------------------------
locus = "chr1:10000-10100"
resp = igv.goto(locus)
print(f"Go to locus response: {resp}")

# -----------------------------
# 5. Snapshot
# -----------------------------
snapshot_dir = tmp_dir / "snapshots"
snapshot_dir.mkdir(exist_ok=True)
resp = igv.snapshot_directory(str(snapshot_dir))
print(f"Snapshot directory response: {resp}")

# -----------------------------
# 6. Create a small network and visualize it
# -----------------------------
network_data = {
    "nodes": [
        {"id": "GeneA", "group": 1},
        {"id": "GeneB", "group": 1},
        {"id": "GeneC", "group": 2},
    ],
    "edges": [
        {"source": "GeneA", "target": "GeneB", "weight": 0.5},
        {"source": "GeneB", "target": "GeneC", "weight": 0.8},
    ]
}

network_viz = NetworkViz()
output_png = tmp_dir / "network_test.png"
network_viz.plot_network(network_data, save_path=str(output_png))
print(f"Network saved to: {output_png.resolve()}")

