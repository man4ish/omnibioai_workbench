import os
import pickle
import json
import subprocess
from datetime import datetime

# -------------------------
# Pickle save/load utilities
# -------------------------
def save_pickle(obj, path: str, overwrite=False):
    """Save a Python object to a .pkl file."""
    if os.path.exists(path) and not overwrite:
        print(f"{path} already exists. Skipping save.")
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(obj, f)
    print(f"Saved object to {path}")
    return path

def load_pickle(path: str):
    """Load a Python object from a .pkl file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} does not exist.")
    with open(path, "rb") as f:
        obj = pickle.load(f)
    return obj

# -------------------------
# Metadata utilities
# -------------------------
def save_metadata(meta: dict, path: str):
    """Save metadata dictionary to a JSON file alongside the model."""
    with open(path, "w") as f:
        json.dump(meta, f, indent=4)
    print(f"Saved metadata to {path}")

def load_metadata(path: str):
    """Load metadata JSON file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} does not exist.")
    with open(path, "r") as f:
        meta = json.load(f)
    return meta

# -------------------------
# Versioned save
# -------------------------
def save_pickle_versioned(
    obj,
    output_dir="models",
    base_name="model",
    description="",
    framework="",
    hyperparameters=None
):
    """
    Save object with timestamp-based versioning and metadata.
    Returns tuple of (pickle_path, metadata_path)
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pickle_path = os.path.join(output_dir, f"{base_name}_{timestamp}.pkl")
    save_pickle(obj, pickle_path)

    # Save metadata
    meta = {
        "model_name": base_name,
        "timestamp": timestamp,
        "description": description,
        "framework": framework,
        "hyperparameters": hyperparameters or {}
    }
    meta_path = os.path.join(output_dir, f"{base_name}_{timestamp}.json")
    save_metadata(meta, meta_path)

    return pickle_path, meta_path

# -------------------------
# DVC integration helpers
# -------------------------
def dvc_add(path: str):
    """Add a file to DVC tracking."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} does not exist.")
    subprocess.run(["dvc", "add", path], check=True)
    print(f"Added {path} to DVC.")

def dvc_push():
    """Push all DVC-tracked files to remote storage."""
    subprocess.run(["dvc", "push"], check=True)
    print("Pushed DVC-tracked files to remote.")

def save_pickle_with_dvc(
    obj,
    output_dir="models",
    base_name="model",
    description="",
    framework="",
    hyperparameters=None
):
    """
    Save an object with versioning, metadata, and add it to DVC.
    Returns tuple of (pickle_path, metadata_path)
    """
    pickle_path, meta_path = save_pickle_versioned(
        obj, output_dir, base_name, description, framework, hyperparameters
    )
    dvc_add(pickle_path)
    dvc_add(meta_path)
    return pickle_path, meta_path
