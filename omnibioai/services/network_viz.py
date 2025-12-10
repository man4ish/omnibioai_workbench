import networkx as nx
import matplotlib.pyplot as plt
import os
from .logger_service import logger

class NetworkVisualizer:
    def __init__(self, output_dir="data/reports"):
        os.makedirs(output_dir, exist_ok=True)
        self.graph = nx.Graph()
        self.output_dir = output_dir

    def build_graph(self, nodes_edges):
        for node in nodes_edges.get("nodes", []):
            self.graph.add_node(node["id"], label=node.get("label", ""))
        for edge in nodes_edges.get("edges", []):
            self.graph.add_edge(edge["source"], edge["target"])
        logger.info(f"Graph built with {self.graph.number_of_nodes()} nodes")

    def save_graph(self, filename="network.png"):
        plt.figure(figsize=(8,6))
        pos = nx.spring_layout(self.graph)
        labels = nx.get_node_attributes(self.graph,'label')
        nx.draw(self.graph, pos, with_labels=True, node_color='skyblue', node_size=1500)
        nx.draw_networkx_labels(self.graph, pos, labels)
        path = os.path.join(self.output_dir, filename)
        plt.savefig(path)
        plt.close()
        logger.info(f"Network graph saved at {path}")
        return path

