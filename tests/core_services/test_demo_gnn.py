import numpy as np
from omnibioai.services.ml_service.gnn_service import GraphData, GNNService

# Simple graph
nodes = 3
edges = [(0,1),(1,0),(1,2),(2,1)]
features = np.eye(3)
graph = GraphData(nodes=nodes, edges=edges, node_features=features)

gnn = GNNService(model_type="gcn", embedding_dim=4, hidden_dim=8, device="cpu")
embeddings = gnn.fit_transform(graph, epochs=20)

print("Node embeddings:")
for i, z in enumerate(embeddings):
    print(f"Node {i}: {z}")

explanation = gnn.explain(embeddings, top_k=1)
print("\nNearest neighbors (top 1):")
for node, neighbors in explanation["similar_nodes"].items():
    print(f"Node {node} -> {neighbors}")
