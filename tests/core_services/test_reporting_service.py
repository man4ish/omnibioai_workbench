import os
import pytest
import pandas as pd
from unittest.mock import AsyncMock, MagicMock
from omnibioai.services.reporting_service import ReportingService

@pytest.fixture
def sample_dataframe():
    return pd.DataFrame({
        "Gene": ["BRCA1", "TP53", "EGFR"],
        "Expression": [12.5, 7.8, 15.2]
    })

@pytest.fixture
def reporting_service(tmp_path):
    # Mock LLMService
    llm_mock = MagicMock()
    llm_mock.generate_async = AsyncMock(return_value="This is a mock summary")
    # Mock NetworkViz and IGVService
    service = ReportingService(report_dir=tmp_path, llm=llm_mock)
    service.network_viz.plot_network = MagicMock()
    service.igv_service.capture_snapshot = MagicMock()
    return service

def test_save_json(reporting_service, tmp_path):
    data = {"key": "value"}
    path = reporting_service.save_json(data, "test.json")
    assert os.path.exists(path)

def test_save_table(reporting_service, sample_dataframe):
    path = reporting_service.save_table(sample_dataframe, "test.csv")
    assert os.path.exists(path)

def test_create_chart(reporting_service, sample_dataframe):
    path = reporting_service.create_chart(sample_dataframe, x_col="Gene", y_col="Expression", filename="test.png")
    assert os.path.exists(path)

@pytest.mark.asyncio
async def test_generate_llm_summary(reporting_service, sample_dataframe):
    summary = await reporting_service.generate_llm_summary(sample_dataframe, context="test context")
    assert summary == "This is a mock summary"

def test_generate_network_figure(reporting_service):
    network_data = {"nodes": [], "edges": []}
    path = reporting_service.generate_network_figure(network_data, filename="network.png")
    assert os.path.exists(path)
    reporting_service.network_viz.plot_network.assert_called_once_with(network_data, save_path=path)

def test_generate_igv_snapshot(reporting_service):
    igv_session = MagicMock()
    path = reporting_service.generate_igv_snapshot(igv_session, filename="igv.png")
    assert os.path.exists(path)
    reporting_service.igv_service.capture_snapshot.assert_called_once_with(igv_session, path)

def test_save_pdf(reporting_service, sample_dataframe):
    pdf_path = reporting_service.save_pdf(
        sample_dataframe,
        filename="report.pdf",
        chart_cols=("Gene", "Expression"),
        agentic_suggestions=["Step1", "Step2"],
        extra_figures=[],
        llm_summary="Mock summary",
        network_figures=[],
        igv_snapshots=[]
    )
    assert os.path.exists(pdf_path)

