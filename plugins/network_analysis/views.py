import csv
import io
import networkx as nx
from django.shortcuts import render
from django.http import HttpResponse
from .forms import UploadNetworkForm

def upload_network(request):
    context = {"form": UploadNetworkForm()}

    if request.method == "POST":
        form = UploadNetworkForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES["file"]
            decoded = uploaded_file.read().decode("utf-8")
            delimiter = '\t' if ".tsv" in uploaded_file.name else ','
            reader = csv.DictReader(io.StringIO(decoded), delimiter=delimiter)

            # Normalize header names
            normalized_headers = {h.lower(): h for h in reader.fieldnames}
            if "source" not in normalized_headers or "target" not in normalized_headers:
                return HttpResponse("Input file must have 'source' and 'target' columns (case-insensitive).", status=400)

            source_key = normalized_headers["source"]
            target_key = normalized_headers["target"]

            G = nx.Graph()
            for row in reader:
                source = row[source_key].strip()
                target = row[target_key].strip()
                if source and target:
                    G.add_edge(source, target)

            nodes = [{"data": {"id": n}} for n in G.nodes()]
            edges = [{"data": {"source": u, "target": v}} for u, v in G.edges()]

            context["elements"] = nodes + edges
            return render(request, "network_analysis/view_network.html", context)

    return render(request, "network_analysis/upload.html", context)
