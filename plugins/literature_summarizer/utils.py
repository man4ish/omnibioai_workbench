from Bio import Entrez
import requests
import json
import os
Entrez.email = "mandecent.gupta@gmail.com"  # Required by NCBI

def fetch_pubmed_abstracts(query, max_results=5):
    """Fetch abstracts from PubMed for a given query."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record.get("IdList", [])
    if not ids:
        return ""
    
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
    abstracts = handle.read()
    handle.close()
    return abstracts



# OLLAMA_ENDPOINT = os.getenv("OLLAMA_ENDPOINT", "http://ollama:11434/api/generate")
# raw_summary = summarize_with_llm(abstracts, endpoint_url=OLLAMA_ENDPOINT)

def summarize_with_llm(text, model="deepseek-r1", temperature=0.7, endpoint_url="http://localhost:11434/api/generate"):
    """Send biomedical text to an LLM API for summarization and stream the response."""
    response = requests.post(
        endpoint_url,
        json={
            "prompt": f"Summarize this biomedical text:\n\n{text}",
            "model": model,
            "temperature": temperature
        },
        stream=True  # Enable streaming for incremental output
    )
    response.raise_for_status()
    
    full_response = ""
    # Stream and accumulate the response line-by-line
    for line in response.iter_lines():
        if line:
            try:
                data = json.loads(line)
                full_response += data.get("response", "")
                if data.get("done"):
                    break
            except json.JSONDecodeError:
                # In case of incomplete line or parse error, skip or handle appropriately
                continue

    print("Full response:", full_response)

    # Attempt to parse final JSON response if available
    try:
        final_data = json.loads(full_response)
    except json.JSONDecodeError:
        # If full response isn't a valid JSON, fallback to dictionary with response text
        final_data = {"response": full_response}

    print("Parsed JSON data:", final_data)  # DEBUG LINE

    if "response" in final_data:
        return final_data["response"]
    else:
        raise ValueError(f"Expected 'response' key not found in API output: {final_data}")

# Optional helper if needed for parsing raw JSON from text
def extract_first_json(text):
    decoder = json.JSONDecoder()
    text = text.lstrip()
    obj, idx = decoder.raw_decode(text)
    return obj, text[idx:].lstrip()
