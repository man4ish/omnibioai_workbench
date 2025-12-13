"""
General-purpose GNN service.

- Input: GraphData (nodes, edges, features)
- Output: node embeddings (and optional artifacts)

Design goals:
- Domain-agnostic
- Backend-pluggable (PyTorch Geometric preferred, graceful fallback)
- Stable public API
"""

from dataclasses import dataclass
from typing import Any, Optional, Dict

import numpy as np

# Optional imports (kept soft to avoid hard dependency failures)
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.data import Data as PyGData
    from torch_geometric.nn import GCNConv, GATConv, SAGEConv
    _HAS_PYG = True
except Exception:
    _HAS_PYG = False


@dataclass
class GraphData:
    """Framework-agnostic graph container."""
    nodes: Any                     # node ids or count
    edges: Any                     # list/array of (src, dst)
    node_features: Optional[Any] = None  # shape: [num_nodes, num_features]
    edge_features: Optional[Any] = None
    graph_features: Optional[Any] = None
    labels: Optional[Any] = None


class _GCN(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int, out_dim: int):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, out_dim)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        return x


class _GAT(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int, out_dim: int, heads: int = 4):
        super().__init__()
        self.conv1 = GATConv(in_dim, hidden_dim, heads=heads)
        self.conv2 = GATConv(hidden_dim * heads, out_dim, heads=1)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.elu(x)
        x = self.conv2(x, edge_index)
        return x


class _GraphSAGE(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int, out_dim: int):
        super().__init__()
        self.conv1 = SAGEConv(in_dim, hidden_dim)
        self.conv2 = SAGEConv(hidden_dim, out_dim)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        return x


class GNNService:
    """
    Public service API.

    Example:
        gnn = GNNService(model_type="gcn", embedding_dim=64)
        Z = gnn.fit_transform(graph_data)
    """

    def __init__(
        self,
        model_type: str = "gcn",      # gcn | gat | graphsage
        embedding_dim: int = 64,
        hidden_dim: int = 128,
        device: Optional[str] = None,
    ):
        if not _HAS_PYG:
            raise ImportError("PyTorch Geometric is required for GNNService")

        self.model_type = model_type.lower()
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")

        self.model: Optional[nn.Module] = None

    # --------- internal helpers ---------

    def _to_pyg(self, graph: GraphData) -> PyGData:
        # Nodes
        if isinstance(graph.nodes, int):
            num_nodes = graph.nodes
        else:
            num_nodes = len(graph.nodes)

        # Node features
        if graph.node_features is None:
            x = torch.eye(num_nodes, dtype=torch.float)
        else:
            x = torch.tensor(graph.node_features, dtype=torch.float)

        # Edges: shape [2, num_edges]
        edge_index = torch.tensor(graph.edges, dtype=torch.long).t().contiguous()

        return PyGData(x=x, edge_index=edge_index)

    def _build_model(self, in_dim: int) -> nn.Module:
        if self.model_type == "gcn":
            return _GCN(in_dim, self.hidden_dim, self.embedding_dim)
        if self.model_type == "gat":
            return _GAT(in_dim, self.hidden_dim, self.embedding_dim)
        if self.model_type in {"graphsage", "sage"}:
            return _GraphSAGE(in_dim, self.hidden_dim, self.embedding_dim)
        raise ValueError(f"Unsupported model_type: {self.model_type}")

    # --------- public API ---------

    def fit_transform(self, graph: GraphData, epochs: int = 50, lr: float = 1e-3) -> np.ndarray:
        """
        Train a GNN (self-supervised / representation learning) and return node embeddings.

        Returns:
            np.ndarray: [num_nodes, embedding_dim]
        """
        pyg_data = self._to_pyg(graph).to(self.device)

        self.model = self._build_model(pyg_data.x.shape[1]).to(self.device)
        optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)

        # Simple unsupervised objective: L2 regularization on embeddings
        # (Placeholder; can be swapped with contrastive / link-prediction losses)
        self.model.train()
        for _ in range(epochs):
            optimizer.zero_grad()
            z = self.model(pyg_data.x, pyg_data.edge_index)
            loss = (z ** 2).mean()
            loss.backward()
            optimizer.step()

        self.model.eval()
        with torch.no_grad():
            embeddings = self.model(pyg_data.x, pyg_data.edge_index)

        return embeddings.cpu().numpy()

    def transform(self, graph: GraphData) -> np.ndarray:
        """Generate node embeddings using a trained model."""
        if self.model is None:
            raise RuntimeError("Model not trained. Call fit_transform first.")

        pyg_data = self._to_pyg(graph).to(self.device)
        self.model.eval()
        with torch.no_grad():
            embeddings = self.model(pyg_data.x, pyg_data.edge_index)
        return embeddings.cpu().numpy()

    def explain(self, embeddings: np.ndarray, top_k: int = 5) -> Dict[str, Any]:
        """
        Lightweight interpretation helper based on embedding similarity.
        """
        # Cosine similarity matrix
        norm = np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-9
        sim = (embeddings @ embeddings.T) / (norm @ norm.T)

        neighbors = {
            i: np.argsort(sim[i])[::-1][1 : top_k + 1].tolist()
            for i in range(sim.shape[0])
        }

        return {
            "similar_nodes": neighbors,
            "embedding_dim": embeddings.shape[1],
        }
