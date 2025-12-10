import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .logger_service import logger
from omnibioai.core.config import REPORT_DIR

class ReportingService:
    def __init__(self, report_dir=REPORT_DIR):
        os.makedirs(report_dir, exist_ok=True)
        self.report_dir = report_dir

    # ----------------------------
    # Save JSON data
    # ----------------------------
    def save_json(self, data, filename="report.json"):
        path = os.path.join(self.report_dir, filename)
        with open(path, "w") as f:
            json.dump(data, f, indent=4)
        logger.info(f"JSON report saved at {path}")
        return path

    # ----------------------------
    # Save Table as CSV
    # ----------------------------
    def save_table(self, dataframe: pd.DataFrame, filename="report.csv"):
        path = os.path.join(self.report_dir, filename)
        dataframe.to_csv(path, index=False)
        logger.info(f"CSV table saved at {path}")
        return path

    # ----------------------------
    # Create chart from dataframe
    # ----------------------------
    def create_chart(self, dataframe: pd.DataFrame, x_col, y_col, chart_type="bar", filename="chart.png"):
        plt.figure(figsize=(8,6))
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

    # ----------------------------
    # Save both table and chart to PDF
    # ----------------------------
    def save_pdf(self, dataframe: pd.DataFrame, filename="report.pdf", chart_cols=None):
        path = os.path.join(self.report_dir, filename)
        with PdfPages(path) as pdf:
            # Table page
            fig, ax = plt.subplots(figsize=(12, len(dataframe)*0.3 + 1))
            ax.axis('tight')
            ax.axis('off')
            table = ax.table(cellText=dataframe.values, colLabels=dataframe.columns, loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            pdf.savefig(fig)
            plt.close(fig)

            # Chart page
            if chart_cols:
                fig, ax = plt.subplots(figsize=(8,6))
                x_col, y_col = chart_cols
                ax.bar(dataframe[x_col], dataframe[y_col])
                ax.set_xlabel(x_col)
                ax.set_ylabel(y_col)
                ax.set_title(f"{y_col} vs {x_col}")
                pdf.savefig(fig)
                plt.close(fig)
        logger.info(f"PDF report saved at {path}")
        return path
