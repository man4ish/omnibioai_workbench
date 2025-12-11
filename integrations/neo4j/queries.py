# Query: genes associated with a specific disease
GET_DISEASE_GENES = """
MATCH (g:Gene)-[:ASSOCIATED_WITH]->(d:Disease {name: $disease})
RETURN g.symbol AS gene, g.name AS name
"""

# Query: network expansion around a gene
EXPAND_GENE_NETWORK = """
MATCH (g:Gene {symbol: $gene})-[r:INTERACTS_WITH]-(n)
RETURN g.symbol AS source, n.symbol AS target, type(r) AS relation
"""

# Query: pathway → gene membership
GET_PATHWAY_GENES = """
MATCH (p:Pathway {name: $pathway})-[:HAS_GENE]->(g:Gene)
RETURN p.name AS pathway, g.symbol AS gene
"""

# Query: co-occurrence between two genes (literature)
GENE_CO_OCCURRENCE = """
MATCH (g1:Gene {symbol: $gene1})-[r:CO_OCCURS]->(g2:Gene {symbol: $gene2})
RETURN r.count AS cooccur_count
"""
