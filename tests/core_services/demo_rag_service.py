# tests/core_services/demo_rag_service.py

from omnibioai.services.rag_service import RAGService

def main():
    service = RAGService(base_output_dir="tests/output_demo")
    query_name = "Alzheimer_CaseStudy"
    query = "Alzheimer Disease AND therapy"

    # Step 1: Download abstracts
    service.download_data(query=query, query_name=query_name, retmax=500, retstart=0)

    # Step 2: Generate embeddings & build FAISS index
    service.generate_embeddings(query_name=query_name)

    # Step 3: Run RAG pipeline
    result = service.run_rag_query(query=query, query_name=query_name, top_k=10, structured=True)

    print("\n=== SUMMARY ===")
    print(result["summary"])
    print("\n=== PMIDs ===")
    print(result["top_pmids"])
    print("\n=== STRUCTURED RESULTS ===")
    print(result["structured_results"])

if __name__ == "__main__":
    main()
