import networkx as nx

def get_sample_network():
    G = nx.karate_club_graph()
    edges = [{'source': u, 'target': v} for u, v in G.edges()]
    nodes = [{'id': n, 'label': str(n)} for n in G.nodes()]
    return {'nodes': nodes, 'edges': edges}

def detect_modules():
    G = nx.karate_club_graph()
    from networkx.algorithms.community import greedy_modularity_communities
    communities = list(greedy_modularity_communities(G))
    module_data = {str(i): list(c) for i, c in enumerate(communities)}
    return {'modules': module_data}