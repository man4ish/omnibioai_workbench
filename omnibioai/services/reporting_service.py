"""
Module: reporting_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the ReportingService class for generating multi-format reports in OmniBioAI.
    Supports JSON and CSV outputs, chart creation, network visualizations, IGV snapshots,
    LLM-generated summaries, and consolidated PDF reports. Designed for scientific reporting
    and data exploration in bioinformatics workflows.

Usage:
    from omnibioai.services.reporting_service import ReportingService
    import pandas as pd
    import asyncio

    df = pd.DataFrame({"Gene": ["GeneA", "GeneB"], "Expression": [10, 15]})
    reporter = ReportingService()

    # Save JSON
    json_path = reporter.save_json({"summary": "Test report"})

    # Save table CSV
    csv_path = reporter.save_table(df)

    # Create chart
    chart_path = reporter.create_chart(df, x_col="Gene", y_col="Expression", chart_type="bar")

    # Generate LLM summary (async)
    llm_summary = asyncio.run(reporter.generate_llm_summary(df, context="Gene expression data"))

    # Save consolidated PDF report
    pdf_path = reporter.save_pdf(
        df,
        filename="final_report.pdf",
        chart_cols=("Gene", "Expression"),
        llm_summary=llm_summary
    )

Classes:
    - ReportingService:
        Service for generating scientific reports in multiple formats. Integrates:
            * LLMService for summaries.
            * NetworkViz for network figures.
            * IGVService for genome browser snapshots.

        Methods:
            * __init__(report_dir=REPORT_DIR, llm: LLMService = None):
                Initializes the reporting service and required components.
            * save_json(data: dict, filename="report.json") -> str:
                Saves a dictionary as a JSON file.
            * save_table(dataframe: pd.DataFrame, filename="report.csv") -> str:
                Saves a DataFrame as a CSV file.
            * create_chart(dataframe: pd.DataFrame, x_col, y_col, chart_type="bar", filename="chart.png") -> str:
                Generates and saves a chart from a DataFrame.
            * generate_llm_summary(dataframe: pd.DataFrame, context: str = "") -> str:
                Generates an asynchronous LLM-based summary for the dataset.
            * generate_network_figure(network_data, filename="network.png") -> str:
                Generates a network visualization and saves it as an image.
            * generate_igv_snapshot(igv_session, filename="igv_snapshot.png") -> str:
                Captures and saves a genome browser snapshot.
            * save_pdf(dataframe: pd.DataFrame, filename="report.pdf", chart_cols: Tuple[str, str] = None, agentic_suggestions: List[str] = None, extra_figures: List[str] = None, llm_summary: str = None, network_figures: List[str] = None, igv_snapshots: List[str] = None) -> str:
                Creates a multi-page PDF report combining tables, charts, figures, LLM summaries, and suggestions.

Dependencies:
    - os, json, pandas, matplotlib.pyplot, matplotlib.backends.backend_pdf.PdfPages
    - LLMService, NetworkViz, IGVService
    - omnibioai.core.config.REPORT_DIR
    - omnibioai.services.logger_service.logger
"""


import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from typing import List, Dict, Tuple
from .logger_service import logger
from .llm_service import LLMService
from .network_viz import NetworkVisualizer
from .igv_service import IGVService
from omnibioai.core.config import REPORT_DIR


class ReportingService:
    def __init__(self, report_dir=REPORT_DIR, llm: LLMService = None):
        os.makedirs(report_dir, exist_ok=True)
        self.report_dir = report_dir
        self.llm = llm or LLMService()
        self.network_viz = NetworkVisualizer()
        self.igv_service = IGVService()

    def save_json(self, data: Dict, filename="report.json") -> str:
        path = os.path.join(self.report_dir, filename)
        with open(path, "w") as f:
            json.dump(data, f, indent=4)
        logger.info(f"JSON report saved at {path}")
        return path

    def save_table(self, dataframe: pd.DataFrame, filename="report.csv") -> str:
        path = os.path.join(self.report_dir, filename)
        dataframe.to_csv(path, index=False)
        logger.info(f"CSV table saved at {path}")
        return path

    def create_chart(self, dataframe: pd.DataFrame, x_col, y_col, chart_type="bar", filename="chart.png") -> str:
        plt.figure(figsize=(8, 6))
        if chart_type == "bar":
            plt.bar(dataframe[x_col], dataframe[y_col])
        elif chart_type == "line":
            plt.plot(dataframe[x_col], dataframe[y_col], marker='o')
        elif chart_type == "scatter":
            plt.scatter(dataframe[x_col], dataframe[y_col])
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title(f"{y_col} vs {x_col}")
        path = os.path.join(self.report_dir, filename)
        plt.savefig(path)
        plt.close()
        logger.info(f"Chart saved at {path}")
        return path

    async def generate_llm_summary(self, dataframe: pd.DataFrame, context: str = "") -> str:
        input_text = (
            f"Summarize the following dataset for a scientific report:\n"
            f"{dataframe.head(10).to_dict()}\nContext: {context}"
        )
        summary = await self.llm.generate_async(input_text)
        logger.info("LLM summary generated")
        return summary

    def generate_network_figure(self, network_data, filename="network.png") -> str:
        path = os.path.join(self.report_dir, filename)
        self.network_viz.plot_network(network_data, save_path=path)
        logger.info(f"Network visualization saved at {path}")
        return path

    def generate_igv_snapshot(self, igv_session, filename="igv_snapshot.png") -> str:
        path = os.path.join(self.report_dir, filename)
        self.igv_service.capture_snapshot(igv_session, path)
        logger.info(f"IGV snapshot saved at {path}")
        return path

    def save_pdf(
        self,
        dataframe: pd.DataFrame,
        filename="report.pdf",
        chart_cols: Tuple[str, str] = None,
        agentic_suggestions: List[str] = None,
        extra_figures: List[str] = None,
        llm_summary: str = None,
        network_figures: List[str] = None,
        igv_snapshots: List[str] = None
    ) -> str:
        path = os.path.join(self.report_dir, filename)
        with PdfPages(path) as pdf:
            # Table
            fig, ax = plt.subplots(figsize=(12, len(dataframe)*0.3 + 1))
            ax.axis("tight")
            ax.axis("off")
            table = ax.table(cellText=dataframe.values, colLabels=dataframe.columns, loc="center")
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            pdf.savefig(fig)
            plt.close(fig)

            # Chart
            if chart_cols:
                x_col, y_col = chart_cols
                fig, ax = plt.subplots(figsize=(8, 6))
                ax.bar(dataframe[x_col], dataframe[y_col])
                ax.set_xlabel(x_col)
                ax.set_ylabel(y_col)
                ax.set_title(f"{y_col} vs {x_col}")
                pdf.savefig(fig)
                plt.close(fig)

            # Agentic AI suggestions
            if agentic_suggestions:
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.axis("off")
                text = "\n".join(f"- {s}" for s in agentic_suggestions)
                ax.text(0.05, 0.95, f"Agentic AI Suggestions:\n{text}", verticalalignment="top", fontsize=10)
                pdf.savefig(fig)
                plt.close(fig)

            # Extra figures
            if extra_figures:
                for fig_path in extra_figures:
                    img = plt.imread(fig_path)
                    fig, ax = plt.subplots(figsize=(8, 6))
                    ax.imshow(img)
                    ax.axis("off")
                    pdf.savefig(fig)
                    plt.close(fig)

            # Network figures
            if network_figures:
                for net_path in network_figures:
                    img = plt.imread(net_path)
                    fig, ax = plt.subplots(figsize=(8, 6))
                    ax.imshow(img)
                    ax.axis("off")
                    pdf.savefig(fig)
                    plt.close(fig)

            # IGV snapshots
            if igv_snapshots:
                for igv_path in igv_snapshots:
                    img = plt.imread(igv_path)
                    fig, ax = plt.subplots(figsize=(8, 6))
                    ax.imshow(img)
                    ax.axis("off")
                    pdf.savefig(fig)
                    plt.close(fig)

            # LLM summary
            if llm_summary:
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.axis("off")
                ax.text(0.05, 0.95, f"LLM Summary:\n\n{llm_summary}", verticalalignment="top", fontsize=10)
                pdf.savefig(fig)
                plt.close(fig)

        logger.info(f"PDF report saved at {path}")
        return path
