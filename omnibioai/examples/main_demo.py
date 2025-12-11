# main_demo.py

import pandas as pd
from omnibioai.services.upload_service import UploadService
#from omnibioai.services.rag_service import RAGServiceCore
from omnibioai.services.llm_service import LLMService
from omnibioai.services.reporting_service import ReportingService
from omnibioai.services.network_viz import NetworkVisualizer
from omnibioai.services.logger_service import logger

# ----------------------------
# 1. Upload Example File
# ----------------------------
upload_service = UploadService()
file_path = upload_service.save_file(
    "example.vcf", 
    b"##fileformat=VCF\n#CHROM\tPOS\tID\nchr1\t123\t.\n"
)
logger.info(f"Uploaded file path: {file_path}")

# ----------------------------
# 2. Run RAG Query
# ----------------------------
#rag_service = RAGServiceCore()
#query = "Genes linked to Parkinson's disease"
#rag_results = rag_service.query(query)
#logger.info(f"RAG query results: {rag_results}")

# ----------------------------
# 3. Generate LLM Summary
# ----------------------------
llm_service = LLMService()
# Combine abstracts/citations as prompt (mock example)
# prompt_text = "Summarize the following citations:\n" + "\n".join(rag_results.get("citations", []))
#summary_text = llm_service.prompt(prompt_text)
#logger.info(f"LLM Summary: {summary_text}")

# ----------------------------
# 4. Create Reporting (Table + Chart + PDF)
# ----------------------------
report_service = ReportingService()

# Prepare DataFrame for reporting
report_data = []
#for i, pmid in enumerate(rag_results.get("citations", [])):
#    report_data.append({"PMID": pmid, "Rank": i+1})
#df = pd.DataFrame(report_data)

# Save JSON
#report_service.save_json(rag_results, filename="rag_results.json")

# Save table CSV
#report_service.save_table(df, filename="rag_results.csv")

# Save chart PNG
report_service.create_chart(df, x_col="PMID", y_col="Rank", chart_type="bar", filename="rag_chart.png")

# Save PDF with table + chart
report_service.save_pdf(df, filename="rag_report.pdf", chart_cols=("PMID","Rank"))

# ----------------------------
# 5. Generate Network Visualization
# ----------------------------
network_service = NetworkVisualizer()
# Mock nodes and edges from RAG results
nodes_edges = {
    "nodes": [{"id": pmid, "label": f"PMID:{pmid}"} for pmid in rag_results.get("citations", [])],
    "edges": [{"source": rag_results.get("citations", [])[i], "target": rag_results.get("citations", [])[i-1]} 
              for i in range(1, len(rag_results.get("citations", [])))]
}
network_service.build_graph(nodes_edges)
network_service.save_graph(filename="rag_network.png")

logger.info("OmniBioAI Core Services Demo Completed Successfully!")

