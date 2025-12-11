from omnibioai.integrations.neo4j.client import Neo4jClient
from omnibioai.integrations.neo4j.queries import EXPAND_GENE_NETWORK
from omnibioai.integrations.neo4j.graph_adapter import GraphAdapter

neo = Neo4jClient(uri, user, password)

records = neo.run_query(EXPAND_GENE_NETWORK, {"gene": "BRCA1"})
graph = GraphAdapter.to_networkx(records)

# Now plugins can visualize or analyze the graph
