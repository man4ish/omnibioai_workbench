# omnibioai/examples/run_agentic_pipeline_example.py

import asyncio
from omnibioai.services.agentic_ai_service import AgenticAIService

async def main():
    # Initialize the Agentic AI Service
    agent_service = AgenticAIService()

    # Example dataset metadata
    dataset_metadata = {
        "name": "20241210_sample_variants",
        "type": "variant",           # 'variant', 'single_cell', or 'multi_omics'
        "task_type": "mutation_analysis",
        "samples": 10
    }

    # List of already run plugin outputs
    plugin_outputs = ["variant_annotation"]  # example: you already ran variant annotation

    # Optional query for gene/pathway/literature suggestions
    query = "BRCA1 related pathways in breast cancer"

    # Run the Agentic AI pipeline
    suggestions = await agent_service.run_agentic_pipeline(
        dataset_metadata=dataset_metadata,
        plugin_outputs=plugin_outputs,
        query=query
    )

    # Print the structured suggestions
    print("=== Agentic AI Suggestions ===")
    for key, value in suggestions.items():
        print(f"{key}: {value}")

# Run asynchronously
if __name__ == "__main__":
    asyncio.run(main())
