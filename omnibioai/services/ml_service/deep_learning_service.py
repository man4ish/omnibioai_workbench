"""
Module: deep_learning_service
Author: Manish Kumar
Version: 1.3
Date: 2025-12-12

Description:
    DeepLearningService class with batch training support, multi-epoch training,
    automatic device detection, and model saving/loading for CNN, LSTM, and Transformers.
"""

from typing import Any, Dict, Optional
import os
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim

class DeepLearningService:
    """
    Service for deep learning models: CNN, LSTM, Transformers.
    """

    def __init__(self, device: Optional[str] = None, model_dir: str = "models"):
        # Automatic device detection
        if device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        elif device not in ["cpu", "cuda"]:
            raise ValueError("Device must be 'cpu' or 'cuda'")
        else:
            self.device = device

        self.model_dir = model_dir
        os.makedirs(self.model_dir, exist_ok=True)
        print(f"[DeepLearningService] Using device: {self.device}")

    # -------------------------
    # Helper: save model
    # -------------------------
    def save_model(self, model: Any, name: str, version: str = "v1") -> str:
        path = os.path.join(self.model_dir, f"{name}_{version}.pt")
        torch.save(model.state_dict() if hasattr(model, "state_dict") else model, path)
        print(f"[DeepLearningService] Model saved: {path}")
        return path

    # -------------------------
    # Helper: create mini-batches
    # -------------------------
    @staticmethod
    def get_batches(X: torch.Tensor, y: torch.Tensor, batch_size: int):
        n_samples = X.shape[0]
        for i in range(0, n_samples, batch_size):
            yield X[i:i+batch_size], y[i:i+batch_size]

    # -------------------------


    # -------------------------
    # CNN Training
    # -------------------------
    def train_cnn(self, X_train: np.ndarray, y_train: np.ndarray, epochs: int = 5,
                lr: float = 0.001, batch_size: int = 16, save: bool = False, version: str = "v1") -> Dict[str, Any]:

        input_size = X_train.shape[1:]
        model = nn.Sequential(
            nn.Conv2d(input_size[0], 16, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Flatten(),
            nn.Linear(16 * input_size[1] * input_size[2], 10)
        ).to(self.device)

        optimizer = optim.Adam(model.parameters(), lr=lr)
        criterion = nn.CrossEntropyLoss()

        X_tensor = torch.tensor(X_train, dtype=torch.float32).to(self.device)
        y_tensor = torch.tensor(y_train, dtype=torch.long).to(self.device)

        model.train()
        for epoch in range(epochs):
            batches = list(self.get_batches(X_tensor, y_tensor, batch_size))
            for batch_idx, (X_batch, y_batch) in enumerate(batches, start=1):
                optimizer.zero_grad()
                outputs = model(X_batch)
                loss = criterion(outputs, y_batch)
                loss.backward()
                optimizer.step()
                print(f"[CNN][Epoch {epoch+1}/{epochs}][Batch {batch_idx}/{len(batches)}] Loss: {loss.item():.4f}")

        if save:
            self.save_model(model, "cnn", version)

        return {"model": model, "status": "trained"}

    # -------------------------
    # LSTM Training
    # -------------------------
    def train_lstm(self, X_train: np.ndarray, y_train: np.ndarray, hidden_size: int = 32,
                epochs: int = 5, lr: float = 0.001, batch_size: int = 16, save: bool = False, version: str = "v1") -> Dict[str, Any]:

        input_size = X_train.shape[2]
        lstm = nn.LSTM(input_size=input_size, hidden_size=hidden_size, batch_first=True).to(self.device)
        linear = nn.Linear(hidden_size, len(np.unique(y_train))).to(self.device)

        optimizer = optim.Adam(list(lstm.parameters()) + list(linear.parameters()), lr=lr)
        criterion = nn.CrossEntropyLoss()

        X_tensor = torch.tensor(X_train, dtype=torch.float32).to(self.device)
        y_tensor = torch.tensor(y_train, dtype=torch.long).to(self.device)

        lstm.train()
        linear.train()
        for epoch in range(epochs):
            batches = list(self.get_batches(X_tensor, y_tensor, batch_size))
            for batch_idx, (X_batch, y_batch) in enumerate(batches, start=1):
                optimizer.zero_grad()
                out, _ = lstm(X_batch)
                out = linear(out[:, -1, :])
                loss = criterion(out, y_batch)
                loss.backward()
                optimizer.step()
                print(f"[LSTM][Epoch {epoch+1}/{epochs}][Batch {batch_idx}/{len(batches)}] Loss: {loss.item():.4f}")

        if save:
            self.save_model((lstm, linear), "lstm", version)

        return {"model": (lstm, linear), "status": "trained"}


    # -------------------------
    # Transformer Training
    # -------------------------
    def train_transformer(self, X_train, y_train, model_name: str = "distilbert-base-uncased",
                          epochs: int = 1, batch_size: int = 4, save: bool = False, version: str = "v1") -> Dict[str, Any]:
        try:
            from transformers import AutoModelForSequenceClassification, AutoTokenizer
        except ImportError:
            return {"status": "transformers library not installed"}

        model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=len(np.unique(y_train)))
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model.to(self.device)

        labels = torch.tensor(y_train, dtype=torch.long).to(self.device)
        optimizer = optim.Adam(model.parameters(), lr=1e-4)
        criterion = nn.CrossEntropyLoss()

        model.train()
        n_samples = len(X_train)
        for epoch in range(epochs):
            for i in range(0, n_samples, batch_size):
                batch_texts = X_train[i:i+batch_size]
                batch_labels = labels[i:i+batch_size]

                encoding = tokenizer(batch_texts, return_tensors="pt", padding=True, truncation=True)
                for key in encoding:
                    encoding[key] = encoding[key].to(self.device)

                optimizer.zero_grad()
                outputs = model(**encoding)
                loss = criterion(outputs.logits, batch_labels)
                loss.backward()
                optimizer.step()
            print(f"[Transformer] Epoch {epoch+1}, Loss: {loss.item():.4f}")

        if save:
            self.save_model(model, "transformer", version)

        return {"model": model, "status": "trained"}
