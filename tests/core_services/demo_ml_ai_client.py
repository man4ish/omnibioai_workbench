# tests/core_services/demo_ml_ai_client.py

import sys
import os

# -------------------------
# Add project root to sys.path
# -------------------------
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))

import pandas as pd
from omnibioai.services.ml_service.ml_ai_client import MLAIClient

def main():
    # -------------------------
    # Sample dataset
    # -------------------------
    data = {
        "feature1": [1, 2, 3, 4, 5, 6],
        "feature2": [10, 20, 30, 40, 50, 60],
        "text": ["good", "bad", "good", "bad", "good", "bad"],
        "label": [0, 1, 0, 1, 0, 1]
    }
    df = pd.DataFrame(data)

    # -------------------------
    # Initialize MLAIClient
    # -------------------------
    client = MLAIClient(use_mlflow=False, device="cpu")

    # -------------------------
    # Rule-based suggestion
    # -------------------------
    suggestion_rule = client.suggest_model(df, goal="classification", mode="rule")
    print("Rule-based model suggestion:")
    print(suggestion_rule)

    # -------------------------
    # AI-based suggestion (LLM placeholder)
    # -------------------------
    suggestion_ai = client.suggest_model(df, goal="classification", mode="ai")
    print("AI-based model suggestion:")
    print(suggestion_ai)

    # -------------------------
    # Train model (rule-based)
    # -------------------------
    print("\nTraining model...")
    result = client.train_model(df, target_col="label", goal="classification", mode="rule")
    print("\nTraining metrics:")
    print(result["metrics"])

if __name__ == "__main__":
    main()
