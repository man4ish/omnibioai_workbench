import os

# ----------------------------
# Workspace / Base Paths
# ----------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
UPLOAD_DIR = os.path.join(DATA_DIR, "uploads")
RESULTS_DIR = os.path.join(DATA_DIR, "results")
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
# Dataset / Reference Data Settings
# ----------------------------
DATASET_DIR = os.path.join(DATA_DIR, "datasets")
DEFAULT_DATASETS = {
    "gnomAD": {"version": "v3.1", "path": os.path.join(DATASET_DIR, "gnomad")},
    "ClinVar": {"version": "202512", "path": os.path.join(DATASET_DIR, "clinvar")},
    "RefSeq": {"version": "2025-12", "path": os.path.join(DATASET_DIR, "refseq")}
}

# ----------------------------
# Annotation Service
# ----------------------------
ANNOTATION_CACHE_SIZE = int(os.environ.get("ANNOTATION_CACHE_SIZE", 10000))

# ----------------------------
# ML / AI Service
# ----------------------------
ML_MODELS_DIR = os.path.join(DATA_DIR, "ml_models")
DEFAULT_ML_MODELS = {
    "variant_predictor": {"version": "1.0", "file": os.path.join(ML_MODELS_DIR, "variant_predictor.pkl")},
    "gene_expression_model": {"version": "1.0", "file": os.path.join(ML_MODELS_DIR, "gene_expression.pkl")}
}

# ----------------------------
# Visualization Service
# ----------------------------
VISUALIZATION_DIR = os.path.join(DATA_DIR, "visualizations")
PLOT_FORMATS = ["png", "pdf", "html"]

# ----------------------------
# Helper function to ensure directories exist
# ----------------------------
def ensure_directories():
    paths = [
        DATA_DIR, UPLOAD_DIR, RESULTS_DIR, REPORT_DIR, LOG_DIR, STATIC_DIR,
        DATASET_DIR, ML_MODELS_DIR, VISUALIZATION_DIR
    ]
    for path in paths:
        os.makedirs(path, exist_ok=True)
