import io
import pandas as pd
import requests
from django.shortcuts import render
from .forms import GeneAnnotationForm
import plotly.express as px
import plotly.io as pio


def ask_llm(prompt: str, model: str, temperature: float) -> str:
    # Call your local Ollama or LLM server (adjust the URL and payload as needed)
    response = requests.post("http://localhost:11434/api/generate", json={
        "model": model,
        "prompt": prompt,
        "temperature": temperature,
        "stream": False
    })
    response.raise_for_status()
    return response.json().get("response", "")


def generate_summary_and_chart(gene_file, prompt, chart_type, model, temperature):
    contents = gene_file.read()
    df = pd.read_csv(io.BytesIO(contents))

    # Generate basic summary stats
    gene_means = df.mean(numeric_only=True).round(2).to_dict()
    stats_summary = "\n".join([f"{k}: {v}" for k, v in gene_means.items()])

    # Prepare prompt for LLM
    llm_prompt = f"This is gene expression data:\n{df.head().to_csv(index=False)}\nSummary statistics:\n{stats_summary}\nUser prompt:\n{prompt}"

    llm_summary = ask_llm(llm_prompt, model, temperature)

    # Prepare data for plotting
    if "Gene" not in df.columns:
        df.insert(0, "Gene", [f"Gene_{i}" for i in range(len(df))])

    melted = df.melt(id_vars=["Gene"], var_name="Sample", value_name="Expression")

    if chart_type == "box":
        fig = px.box(melted, x="Sample", y="Expression", color="Sample", title="Expression Distribution per Sample")
    elif chart_type == "scatter":
        fig = px.scatter(melted, x="Sample", y="Expression", color="Sample", title="Expression Scatter Plot")
    elif chart_type == "violin":
        fig = px.violin(melted, x="Sample", y="Expression", color="Sample", box=True, points="all", title="Violin Plot")
    elif chart_type == "heatmap":
        pivot_df = df.set_index("Gene")
        fig = px.imshow(pivot_df.T, labels=dict(x="Gene", y="Sample", color="Expression"), title="Gene Expression Heatmap")
    else:
        fig = px.box(melted, x="Sample", y="Expression", color="Sample", title="Default Box Plot")

    chart_html = pio.to_html(fig, include_plotlyjs='cdn')
    return llm_summary, chart_html


def annotate_view(request):
    if request.method == "POST":
        form = GeneAnnotationForm(request.POST, request.FILES)
        if form.is_valid():
            summary, chart_html = generate_summary_and_chart(
                form.cleaned_data['gene_file'],
                form.cleaned_data['prompt'],
                form.cleaned_data['chart_type'],
                form.cleaned_data['model'],
                form.cleaned_data['temperature'],
            )
            return render(request, "gene_annotation/result.html", {
                "summary": summary,
                "chart_html": chart_html,
            })
    else:
        form = GeneAnnotationForm()

    return render(request, "gene_annotation/upload.html", {"form": form})
