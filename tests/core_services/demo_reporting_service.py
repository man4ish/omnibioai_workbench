"""
Demo script for OmniBioAI ReportingService.

Demonstrates:
1. Saving JSON report
2. Saving DataFrame as CSV
3. Creating charts
4. Generating LLM summary (async)
5. Creating network visualization
6. Capturing IGV snapshot
7. Consolidating all components into a PDF report
"""

import pandas as pd
import asyncio
from omnibioai.services.reporting_service import ReportingService
from omnibioai.services.igv_helpers import write_igv_html

def main():
    # Initialize ReportingService
    reporter = ReportingService(report_dir="tests/output_reporting_demo")

    # Example DataFrame
    df = pd.DataFrame({
        "Gene": ["GeneA", "GeneB", "GeneC"],
        "Expression": [10, 15, 7],
        "PValue": [0.01, 0.05, 0.20]
    })

    print("\n=== Step 1: Save JSON ===")
    json_path = reporter.save_json({"summary": "Test report", "dataset_size": len(df)})
    print(f"JSON saved at: {json_path}")

    print("\n=== Step 2: Save CSV table ===")
    csv_path = reporter.save_table(df)
    print(f"CSV saved at: {csv_path}")

    print("\n=== Step 3: Create chart ===")
    chart_path = reporter.create_chart(df, x_col="Gene", y_col="Expression", chart_type="bar")
    print(f"Chart saved at: {chart_path}")

    print("\n=== Step 4: Generate LLM summary (async) ===")
    #llm_summary = asyncio.run(reporter.generate_llm_summary(df, context="Gene expression dataset"))
    #print(f"LLM summary:\n{llm_summary}")

    print("\n=== Step 5: Generate network figure ===")
    # Example network data (replace with actual network object as needed)
    network_data = {
        "nodes": [{"id": "GeneA"}, {"id": "GeneB"}],
        "edges": [{"source": "GeneA", "target": "GeneB"}]
    }
    # Build the graph first
    reporter.network_viz.build_graph(network_data)

    # Save a static PNG
    network_path = reporter.network_viz.save_static_graph("network.png")

    print(f"Network figure saved at: {network_path}")

    print("\n=== Step 6: Generate IGV snapshot ===")
    # Example IGV session
    

    tracks = [
        {"name": "Example BAM", "url": "data/sample.bam", "type": "alignment", "format": "bam"}
    ]

    write_igv_html(
        output_html_path="tests/output_reporting_demo/igv_demo.html",
        genome="hg38",
        tracks=tracks,
        locus="chr1:100000-200000"
    )


    """
    print("\n=== Step 7: Save consolidated PDF report ===")
    pdf_path = reporter.save_pdf(
        df,
        filename="final_report.pdf",
        chart_cols=("Gene", "Expression"),
        llm_summary=llm_summary,
        network_figures=[network_path],
        igv_snapshots=[igv_path]
    )
    print(f"PDF report saved at: {pdf_path}")
    """
 
if __name__ == "__main__":
    main()
