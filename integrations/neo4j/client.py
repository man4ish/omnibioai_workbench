import logging
from neo4j import GraphDatabase
from .exceptions import Neo4jConnectionError, Neo4jQueryError

logger = logging.getLogger(__name__)

class Neo4jClient:
    def __init__(self, uri, user, password):
        try:
            self.driver = GraphDatabase.driver(uri, auth=(user, password))
            logger.info("Connected to Neo4j at %s", uri)
        except Exception as e:
            raise Neo4jConnectionError(f"Could not connect to Neo4j: {e}")

    def close(self):
        if self.driver:
            self.driver.close()

    def run_query(self, cypher, params=None):
        try:
            with self.driver.session() as session:
                result = session.run(cypher, params or {})
                return [record.data() for record in result]
        except Exception as e:
            raise Neo4jQueryError(f"Cypher query failed: {e}")
