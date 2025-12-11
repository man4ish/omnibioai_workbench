import networkx as nx

class GraphAdapter:
    """
    Converts Neo4j Cypher results into graph structures usable by NetworkX.
    """

    @staticmethod
    def to_networkx(edge_records, source_key="source", target_key="target"):
        """
        Convert Neo4j edge records into a NetworkX Graph.
        
        edge_records example:
        [
            {"source": "BRCA1", "target": "BRCA2", "relation": "INTERACTS_WITH"},
            ...
        ]
        """
        G = nx.Graph()
        for rec in edge_records:
            src = rec.get(source_key)
            tgt = rec.get(target_key)
            rel = rec.get("relation", "")

            G.add_node(src)
            G.add_node(tgt)
            G.add_edge(src, tgt, relation=rel)

        return G

    @staticmethod
    def to_pyvis(edge_records):
        """
        Convert Neo4j output to PyVis-compatible dict.
        """
        nodes = set()
        edges = []

        for rec in edge_records:
            src = rec["source"]
            tgt = rec["target"]
            nodes.add(src)
            nodes.add(tgt)

            edges.append({
                "from": src,
                "to": tgt,
                "label": rec.get("relation", "")
            })

        return {
            "nodes": [{"id": n, "label": n} for n in nodes],
            "edges": edges
        }
