import os
import json
import shutil
from omnibioai.services.network_viz import NetworkVisualizer
from omnibioai.services.logger_service import logger

# ----------------------------
# 1. Initialize the visualizer
# ----------------------------
test_output_dir = "data/test_network"
visualizer = NetworkVisualizer(output_dir=test_output_dir)

# ----------------------------
# 2. Define sample graph
# ----------------------------
nodes_edges = {
    "nodes": [
        {"id": "n1", "label": "Node 1"},
        {"id": "n2", "label": "Node 2"},
        {"id": "n3", "label": "Node 3"}
    ],
    "edges": [
        {"source": "n1", "target": "n2"},
        {"source": "n2", "target": "n3"}
    ]
}

# ----------------------------
# 3. Build the graph
# ----------------------------
visualizer.build_graph(nodes_edges)

# ----------------------------
# 4. Save static PNG
# ----------------------------
png_path = visualizer.save_static_graph(filename="test_network.png")
if os.path.exists(png_path):
    logger.info(f"Static graph created successfully: {png_path}")
else:
    logger.error("Failed to create static graph PNG")

# ----------------------------
# 5. Export Cytoscape JSON
# ----------------------------
json_path = visualizer.export_cytoscape_json(filename="test_network.json")
if os.path.exists(json_path):
    logger.info(f"Cytoscape JSON exported successfully: {json_path}")
    with open(json_path, "r") as f:
        data = json.load(f)
        # Handle different NetworkX versions
        edge_key = "links" if "links" in data else "edges"
        assert "nodes" in data and edge_key in data, "Invalid Cytoscape JSON structure"
        logger.info(f"Cytoscape JSON structure verified: {len(data['nodes'])} nodes, {len(data[edge_key])} edges")
else:
    logger.error("Failed to export Cytoscape JSON")

# ----------------------------
# 6. Cleanup
# ----------------------------
shutil.rmtree(test_output_dir)
logger.info(f"Cleaned up test directory: {test_output_dir}")
