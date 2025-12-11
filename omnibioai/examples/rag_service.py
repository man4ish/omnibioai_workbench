from integrations.ragbio import RAGBioClient, GeneDiscoveryPipeline

rag = GeneDiscoveryPipeline(base_dir="data/PubMed")

result = rag.run(
    query="Alzheimer Disease AND therapy",
    query_name="Alzheimer_CaseStudy"
)

print(result["abstracts"][:5])
