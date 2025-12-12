"""
utils_nlp.py
------------

Utility functions for NLP workflows using PyTorch, HuggingFace Transformers, and scikit-learn. 
Includes helpers for reproducibility, device selection, timing, model/tokenizer loading, embedding 
normalization and similarity computations, standard metrics, and miscellaneous utilities.

Main Features
-------------
1. Reproducibility & Device
   - set_seed(seed: int = 42)
       Set random seed for Python, NumPy, and PyTorch for reproducible results.
   - get_device() -> str
       Return 'cuda' if a GPU is available, else 'cpu'.

2. Timer for benchmarking
   - timer(name="task")
       Context manager to measure execution time of a code block.

3. Model & Tokenizer helpers
   - load_hf_model(model_name: str, task: str = None, device: str = None)
       Load a HuggingFace model and tokenizer, or a pipeline for a specified task.

4. Embedding utilities
   - normalize_embeddings(embeddings)
       Normalize embeddings to unit L2 norm.
   - cosine_similarity(a: torch.Tensor, b: torch.Tensor)
       Compute cosine similarity between two embedding matrices.

5. Metrics
   - compute_metrics(y_true, y_pred)
       Compute standard metrics: accuracy and macro F1 score.

6. Miscellaneous Helpers
   - ensure_dir(path: str)
       Create a directory if it does not exist.

Dependencies
------------
- os, random, time, contextlib
- numpy, torch, scikit-learn, transformers

Usage Example
-------------
from utils_nlp import set_seed, get_device, load_hf_model, compute_metrics

set_seed(123)
device = get_device()
model, tokenizer = load_hf_model("bert-base-uncased")
metrics = compute_metrics([0,1], [0,1])
"""


import os
import random
import time
from contextlib import contextmanager

import numpy as np
import torch
from sklearn.preprocessing import normalize
from sklearn.metrics import accuracy_score, f1_score
from transformers import AutoModel, AutoTokenizer, pipeline

# -------------------------
# Reproducibility & Device
# -------------------------
def set_seed(seed: int = 42):
    """Set random seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def get_device():
    """Return 'cuda' if GPU is available, else 'cpu'."""
    return "cuda" if torch.cuda.is_available() else "cpu"


# -------------------------
# Timer for benchmarking
# -------------------------
@contextmanager
def timer(name="task"):
    """Context manager to measure execution time."""
    start = time.time()
    yield
    print(f"[{name}] took {time.time() - start:.2f}s")


# -------------------------
# Model & Tokenizer helpers
# -------------------------
def load_hf_model(model_name: str, task: str = None, device: str = None):
    """
    Load a HuggingFace model and tokenizer or pipeline.
    
    Parameters
    ----------
    model_name : str
        Model identifier.
    task : str, optional
        Pipeline task name (e.g., 'text-classification').
    device : str, optional
        'cuda' or 'cpu'. Defaults to GPU if available.
        
    Returns
    -------
    If task is specified: transformers.Pipeline
    Else: tuple(torch.nn.Module, AutoTokenizer)
    """
    device = device or get_device()
    if task:
        return pipeline(task, model=model_name, tokenizer=model_name, device=0 if device=="cuda" else -1)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to(device)
    return model, tokenizer


# -------------------------
# Embedding utilities
# -------------------------
def normalize_embeddings(embeddings):
    """
    Normalize embeddings to unit L2 norm.
    
    Parameters
    ----------
    embeddings : np.ndarray or torch.Tensor
    
    Returns
    -------
    torch.Tensor
    """
    if isinstance(embeddings, np.ndarray):
        embeddings = torch.tensor(embeddings, dtype=torch.float32)
    return torch.tensor(normalize(embeddings, axis=1, norm='l2'))


def cosine_similarity(a: torch.Tensor, b: torch.Tensor):
    """
    Compute cosine similarity between two embedding matrices.
    
    Parameters
    ----------
    a, b : torch.Tensor [n_samples, dim]
    
    Returns
    -------
    torch.Tensor [n_samples, n_samples]
    """
    a = normalize_embeddings(a)
    b = normalize_embeddings(b)
    return a @ b.T


# -------------------------
# Metrics
# -------------------------
def compute_metrics(y_true, y_pred):
    """
    Compute standard metrics: accuracy and macro F1.
    """
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "f1": f1_score(y_true, y_pred, average="macro")
    }


# -------------------------
# Misc Helpers
# -------------------------
def ensure_dir(path: str):
    """Create a directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)
    return path
