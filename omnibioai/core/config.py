import os

# ----------------------------
# Workspace / Base Paths
# ----------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
UPLOAD_DIR = os.path.join(DATA_DIR, "uploads")
REPORT_DIR = os.path.join(DATA_DIR, "reports")
LOG_DIR = os.path.join(DATA_DIR, "logs")
STATIC_DIR = os.path.join(BASE_DIR, "static")

# ----------------------------
# RAG / Model Settings
# ----------------------------
RAG_EMBEDDING_MODEL = os.environ.get("RAG_EMBEDDING_MODEL", "BioBERT")
RAG_LLM_MODEL = os.environ.get("RAG_LLM_MODEL", "deepseek-r1:latest")
TOP_K = int(os.environ.get("RAG_TOP_K", 5))
USE_GPU = os.environ.get("RAG_USE_GPU", "False").lower() == "true"

# ----------------------------
# Database / Neo4j
# ----------------------------
NEO4J_URI = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.environ.get("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.environ.get("NEO4J_PASSWORD", "password")

# ----------------------------
# File / Upload Settings
# ----------------------------
ALLOWED_FILE_EXTENSIONS = ["vcf", "h5", "csv", "json", "pdf"]

# ----------------------------
# Logging
# ----------------------------
LOG_LEVEL = os.environ.get("LOG_LEVEL", "INFO")

# ----------------------------
# Helper function to ensure directories exist
# ----------------------------
def ensure_directories():
    for path in [DATA_DIR, UPLOAD_DIR, REPORT_DIR, LOG_DIR, STATIC_DIR]:
        os.makedirs(path, exist_ok=True)
