"""
Demo: DeepLearningService with batch training and multi-epoch support
Author: Manish Kumar
Date: 2025-12-12
"""

import sys
import os
import numpy as np

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))

from omnibioai.services.ml_service.deep_learning_service import DeepLearningService

def demo_cnn(dl_service):
    print("=== CNN Demo ===")
    X_train = np.random.randn(64, 3, 32, 32).astype(np.float32)
    y_train = np.random.randint(0, 10, 64)

    result = dl_service.train_cnn(X_train, y_train, epochs=3, batch_size=16, save=True, version="v1")
    print("CNN Status:", result["status"])
    print("CNN Model:", result["model"])
    print()

def demo_lstm(dl_service):
    print("=== LSTM Demo ===")
    X_train = np.random.randn(32, 10, 8).astype(np.float32)
    y_train = np.random.randint(0, 5, 32)

    result = dl_service.train_lstm(X_train, y_train, hidden_size=16, epochs=2, batch_size=8, save=True, version="v1")
    print("LSTM Status:", result["status"])
    print("LSTM Model:", result["model"])
    print()

def demo_transformer(dl_service):
    print("=== Transformer Demo ===")
    X_train = ["This is a test sentence." for _ in range(8)]
    y_train = np.random.randint(0, 2, 8)

    result = dl_service.train_transformer(X_train, y_train, model_name="distilbert-base-uncased",
                                          epochs=1, batch_size=4, save=True, version="v1")
    print("Transformer Status:", result["status"])
    if "model" in result:
        print("Transformer Model:", result["model"])
    print()

if __name__ == "__main__":
    # Initialize service with automatic device detection
    dl_service = DeepLearningService()
    print(f"[Demo] Using device: {dl_service.device}\n")

    demo_cnn(dl_service)
    demo_lstm(dl_service)
    demo_transformer(dl_service)

    print("[Demo] All models trained and saved with batch training and multi-epoch support.")
