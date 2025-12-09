from django.shortcuts import render
from .forms import LiteratureForm
from .utils import fetch_pubmed_abstracts, summarize_with_llm
from collections import defaultdict, Counter
from itertools import combinations
import networkx as nx
import requests
from xml.etree import ElementTree as ET
import json
import re
from django.shortcuts import render
from collections import defaultdict
import requests
from rake_nltk import Rake
import html
import datetime


def extract_keywords(text, min_length=2, max_keywords=10):
    """
    Extracts keywords from text using RAKE.
    
    Parameters:
        text (str): The input abstract or document text.
        min_length (int): Minimum word length to include in keywords.
        max_keywords (int): Maximum number of top keywords to return.
    
    Returns:
        List[str]: A list of keywords.
    """
    rake = Rake()  # Uses NLTK stopwords by default
    rake.extract_keywords_from_text(text)
    ranked_phrases = rake.get_ranked_phrases()
    
    # Filter out too-short keywords and return top-N
    filtered_keywords = [
        phrase.lower() for phrase in ranked_phrases
        if len(phrase.split()) >= min_length
    ][:max_keywords]

    return filtered_keywords

def chunk_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def pubmed_query(keyword, max_results=200):
    """
    Query PubMed for a keyword and return a list of articles with metadata.
    Fetches up to max_results articles in batches of 100.
    """
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    esearch_params = {
        "db": "pubmed",
        "term": keyword,
        "retmax": str(max_results),
        "retmode": "xml",
        "sort": "pub+date"  # Sort by publication date, newest first
    }
    esearch_resp = requests.get(esearch_url, params=esearch_params)
    if esearch_resp.status_code != 200:
        print("ESearch request failed")
        return []

    esearch_root = ET.fromstring(esearch_resp.text)
    pmid_list = [id_elem.text for id_elem in esearch_root.findall(".//Id")]

    if not pmid_list:
        print("No PMIDs found")
        return []

    results = []
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    for pmid_chunk in chunk_list(pmid_list, 100):  # batch size of 100 for EFetch
        efetch_params = {
            "db": "pubmed",
            "id": ",".join(pmid_chunk),
            "retmode": "xml"
        }
        efetch_resp = requests.get(efetch_url, params=efetch_params)
        if efetch_resp.status_code != 200:
            print(f"EFetch failed for PMIDs: {pmid_chunk}")
            continue

        root = ET.fromstring(efetch_resp.text)

        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//PMID")
            title = article.findtext(".//ArticleTitle")

            abstract_text = ""
            abstract_elements = article.findall(".//AbstractText")
            if abstract_elements:
                abstract_text = "".join([elem.text or "" for elem in abstract_elements])

            pub_date = ""
            pubdate_elem = article.find(".//JournalIssue/PubDate")
            if pubdate_elem is not None:
                year = pubdate_elem.findtext("Year")
                month = pubdate_elem.findtext("Month")
                day = pubdate_elem.findtext("Day")

                if year:
                    pub_date = year
                    if month:
                        pub_date += f" {month}"
                    if day:
                        pub_date += f" {day}"
                else:
                    medline_date = pubdate_elem.findtext("MedlineDate")
                    if medline_date:
                        pub_date = medline_date

            results.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract_text,
                "pubDate": pub_date
            })

    return results

def parse_pubtator(pubtator_text):
    """
    Parse PubTator text to extract entity annotations.
    Returns a list of dicts with keys: start, end, type, text.
    """
    annotations = []
    for line in pubtator_text.splitlines():
        parts = line.split('\t') if '\t' in line else line.split(' ')
        if len(parts) >= 5 and parts[0].isdigit():
            # Format: PMID start end type id (PubTator uses spaces for annotations)
            try:
                pmid, start, end, entity_type, entity_id = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
                annotations.append({"start": start, "end": end, "type": entity_type})
            except Exception:
                # Skip malformed lines
                continue
    return annotations


def highlight_entities(text, annotations):
    """
    Given a plain text and list of entity annotations (with start/end),
    return HTML string with entities highlighted.
    """
    if not annotations:
        return html.escape(text)  # Escape HTML chars if no annotations

    # Sort annotations by start index (important for correct slicing)
    annotations = sorted(annotations, key=lambda x: x['start'])

    highlighted = []
    last_idx = 0

    for ann in annotations:
        start, end, ent_type = ann['start'], ann['end'], ann['type']
        # Append text before entity
        highlighted.append(html.escape(text[last_idx:start]))
        # Append highlighted entity span
        entity_text = html.escape(text[start:end])
        highlighted.append(f'<span class="entity {ent_type.lower()}">{entity_text}</span>')
        last_idx = end

    # Append remaining text after last entity
    highlighted.append(html.escape(text[last_idx:]))

    return ''.join(highlighted)

def pubtator_annotate(request):
    query = request.GET.get("query", "")
    annotated = []

    if query:
        results = pubmed_query(query)
        print(f"PubMed query results for '{query}': {results}")  # Debug line
        
        current_year = datetime.datetime.now().year
        annotated = []

        for res in results[:10]:
            pmid = res.get("pmid")
            pub_date = res.get("pubDate", "")
            
            # Skip if no PMID or future publication date
            if not pmid or not pub_date or not pub_date[:4].isdigit() or int(pub_date[:4]) > current_year:
                print(f"Skipping PMID {pmid} due to future or invalid pubDate: {pub_date}")
                continue

            # Call PubTator
            r = requests.get(
                f"https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids={pmid}"
            )
            print(f"PubTator API response for PMID {pmid}: status {r.status_code}")
            print(r.text[:500])  # print first 500 chars of response

            if r.status_code == 200:
                annotations = parse_pubtator(r.text)
                print(f"Parsed {len(annotations)} annotations for PMID {pmid}")
                highlighted = highlight_entities(res["abstract"], annotations)
                annotated.append({
                    "title": res["title"],
                    "text": highlighted
                })

    if not annotated and query:
        # Show a friendly message if no abstracts found
        annotated = None

    return render(request, "literature_summarizer/annotated_text.html", {"results": annotated, "query": query})


def keyword_network(request):
    """
    Generate keyword co-occurrence network graph from abstracts.
    """
    query = request.GET.get("query", "")
    co_occurrences = Counter()

    if query:
        abstracts = [r["abstract"] for r in pubmed_query(query)]
        for abs_text in abstracts:
            keywords = extract_keywords(abs_text)  # Implement or import this (e.g., RAKE, TF-IDF)
            for combo in combinations(set(keywords), 2):
                co_occurrences[tuple(sorted(combo))] += 1

    G = nx.Graph()
    for (k1, k2), weight in co_occurrences.items():
        if weight >= 2:
            G.add_edge(k1, k2, weight=weight)

    nodes = [{"data": {"id": n}} for n in G.nodes()]
    edges = [{"data": {"source": u, "target": v, "weight": d['weight']}} for u, v, d in G.edges(data=True)]

    return render(request, "literature_summarizer/cooccurrence.html", {"elements": nodes + edges})

def extract_year(pub_date_str):
    """
    Extracts a 4-digit year from a date string.
    Returns None if no valid year is found.
    """
    match = re.search(r"\b(19|20)\d{2}\b", pub_date_str)
    return match.group(0) if match else None

def trend_analysis(request):
    keyword = request.GET.get("keyword", "")
    year_freq = defaultdict(int)

    if keyword:
        results = pubmed_query(keyword)
        for result in results:
            year = extract_year(result.get("pubDate", ""))
            if year:
                year_freq[year] += 1

    data = sorted(year_freq.items())  # List of tuples (year, frequency)

    labels = [year for year, freq in data]
    values = [freq for year, freq in data]

    # Pass JSON strings so Chart.js gets proper JS arrays
    context = {
        "keyword": keyword,
        "labels": json.dumps(labels),
        "values": json.dumps(values),
    }

    return render(request, "literature_summarizer/trends.html", context)

    
    return render(request, "literature_summarizer/trends.html", {"data": data, "keyword": keyword})

def clean_llm_output(output: str) -> str:
    """
    Clean LLM output by removing tags like <think> or unwanted metadata.
    """
    lines = output.splitlines()
    cleaned_lines = [line for line in lines if line.strip() and not line.strip().startswith('<')]
    return "\n".join(cleaned_lines).strip()

def summarize_view(request):
    """
    Handle literature summarization form and display summary.
    """
    summary = ""
    if request.method == "POST":
        form = LiteratureForm(request.POST)
        if form.is_valid():
            query = form.cleaned_data["query"]
            num = form.cleaned_data["num_results"]
            abstracts = fetch_pubmed_abstracts(query, max_results=num)
            try:
                import os
                raw_summary = summarize_with_llm(abstracts, endpoint_url="http://localhost:11434/api/generate")
                summary = clean_llm_output(raw_summary)
            except Exception as e:
                summary = f"Error during summarization: {e}"
    else:
        form = LiteratureForm()
    return render(request, "literature_summarizer/summarize.html", {"form": form, "summary": summary})
