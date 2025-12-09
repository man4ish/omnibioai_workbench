# OmniBioAI

A modular, **plugin-based Django bioinformatics workbench** that integrates tools for genomic data visualization, machine learning prediction, variant annotation, pathway enrichment analysis, literature mining, single-cell analysis, and network analysis.

---

## Project Structure

```
omnibioai_project/
├── manage.py                    # Django management script
├── db.sqlite3                   # SQLite DB for development
├── start_app.sh                 # Shell script to run the app
├── data/                        # Input/output data directory
├── plugins/                     # All plugin-based apps
│   ├── home/                    # Landing page plugin
│   ├── igv_viewer/              # Genomic viewer integration (IGV)
│   ├── pipeline_manager/        # Workflow management plugin
│   ├── ml_predictor/            # Machine learning prediction plugin
│   ├── variant_annotation/      # Variant annotation plugin
│   ├── pathway_enrichment/      # Pathway enrichment & GSEA plugin
│   ├── network_analysis/        # Keyword co-occurrence network visualization plugin
│   ├── literature_summarizer/   # Literature mining & summarization plugin
│   ├── single_cell_analysis/    # Single-cell RNA-seq processing & visualization
│   ├── gene_annotation/         # Functional gene and variant annotation plugin
├── omnibioai/                   # Main Django project configuration
│   ├── settings.py
│   ├── urls.py
│   └── ...
├── requirements.txt
└── Dockerfile
```

---

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/omnibioai_project.git
cd omnibioai_project
```

### 2. Set Up Virtual Environment

```bash
python3 -m venv env
source env/bin/activate
```

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```

> If `requirements.txt` is missing, install core dependencies manually:

```bash
pip install django gseapy pandas
```

### 4. Run Migrations

```bash
python manage.py migrate
```

### 5. Start the Server

```bash
./start_app.sh
# or
python manage.py runserver
```

Visit [http://127.0.0.1:8000](http://127.0.0.1:8000) to access the workbench.

---

## Plugins Overview

| Plugin                   | Description                                                                                  |
| ------------------------ | -------------------------------------------------------------------------------------------- |
| `home/`                  | Landing page with general info about OmniBioAI                                               |
| `igv_viewer/`            | Visualize genomic tracks via IGV.js                                                          |
| `pipeline_manager/`      | Launch and monitor bioinformatics pipelines                                                  |
| `ml_predictor/`          | Apply machine learning models to biological datasets                                         |
| `variant_annotation/`    | Annotate VCF or variant files using bioinformatics tools                                     |
| `pathway_enrichment/`    | Run Enrichr or GSEA on gene lists and visualize results                                      |
| `network_analysis/`      | Analyze and visualize keyword co-occurrence networks from literature                         |
| `literature_summarizer/` | Summarize, extract entities, and analyze trends from biomedical literature using LLMs        |
| `single_cell_analysis/`  | Single-cell RNA-seq data processing and visualization                                        |
| `gene_annotation/`       | Functional annotation of genes and variants leveraging LLMs                                  |
| `ollama-server/`         | Backend API for Ollama LLM integration                                                       |
| `bio_navigator/`         | Interactive genomic data query and analysis using Hugging Face-powered AI                    |
| `data_uploader/`         | Upload and process biological datasets for various plugins, including FAISS indexing for RAG |
| `rag_inference/`         | Retrieval-Augmented Generation (RAG) for querying biological documents                       |

---

## Pathway Enrichment Input Formats

* **Enrichr**: `.tsv` file with **one gene per line**
* **GSEA**: `.tsv` file with **two columns**: `gene` and `score`

---

## Utilities

* `start_app.sh`: Launch the Django application
* `data/`: Temporary directory for uploaded and processed files

---

## Docker Usage

### Build Docker Image

```bash
docker build -t omnibioai:latest .
```

### Run Docker Container

```bash
docker run -d -p 8000:8000 omnibioai:latest
```

* `-d`: Run container in background
* `-p 8000:8000`: Map local port 8000 to container port 8000

Visit [http://localhost:8000](http://localhost:8000)

---

## Kubernetes Deployment (Optional)

You can deploy OmniBioAI and Ollama Server using Kubernetes for scalable, production-ready setups.

### Steps:

1. Create namespace:

```bash
kubectl create namespace omnibioai
```

2. Apply manifests:

```bash
kubectl apply -f K8s/omnibioai-deployment.yaml
kubectl apply -f K8s/omnibioai-service.yaml
kubectl apply -f K8s/ollama-deployment.yaml
kubectl apply -f K8s/ollama-service.yaml
```

3. Access services:

* OmniBioAI: `http://localhost:30080` (NodePort example)
* Ollama Server: `http://localhost:31434`

For Minikube:

```bash
minikube service omnibioai-service -n omnibioai
minikube service ollama-service -n omnibioai
```

---

![Screenshot](./screenshots/homepage.png)


