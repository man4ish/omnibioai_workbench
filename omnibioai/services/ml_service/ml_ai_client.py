"""
Module: ml_ai_client
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the MLAIClient class, which acts as a unified interface for
    machine learning, deep learning, NLP, time series, and survival analysis.
    Supports both rule-based and LLM-based AI suggestions for model selection,
    automatic dispatch to specialized services, and optional MLflow tracking.
    Includes utility functions for model evaluation and visualization.

Usage:
    from omnibioai.ml_service.ml_ai_client import MLAIClient
    import pandas as pd

    df = pd.read_csv("dataset.csv")
    client = MLAIClient(use_mlflow=True, device="cuda")

    # Suggest a model
    suggestion = client.suggest_model(df, goal="classification", mode="ai")

    # Train the suggested model
    result = client.train_model(df, target_col="label", goal="classification", mode="ai")
    print(result["metrics"])

Classes:
    - MLAIClient:
        Unified ML client with rule-based and LLM-based model suggestions.

        Attributes:
            deep_learning_service (DeepLearningService): Handles CNN/LSTM/Transformer models.
            time_series_service (TimeSeriesService): Handles time series modeling.
            nlp_service (NLPService): Handles NLP tasks.
            survival_service (SurvivalService): Handles survival analysis.
            mlflow_service (MLFlowTrackingService | None): Optional MLflow logger.
            llm_service (LLMService): Used for AI-based model suggestions.

        Methods:
            suggest_model(df: pd.DataFrame, goal: Optional[str] = None, mode: str = "rule") -> Dict[str, Any]:
                Suggests the best model for a dataset using rules or LLM-based AI.

            train_model(df: pd.DataFrame, target_col: str, goal: Optional[str] = None, mode: str = "rule", **kwargs) -> Dict[str, Any]:
                Suggests and trains a model automatically, dispatching to the appropriate service.
                Optionally logs to MLflow and produces evaluation plots.

Dependencies:
    - pandas: For DataFrame manipulation.
    - DeepLearningService, TimeSeriesService, NLPService, SurvivalService, MLFlowTrackingService
      (from the omnibioai.ml_service package)
    - LLMService: For AI-based model suggestions.
"""

from typing import Optional, Dict, Any
import pandas as pd

from .deep_learning_service import DeepLearningService
from .time_series_service import TimeSeriesService
from .nlp_service import NLPService
from .survival_service import SurvivalService
from .mlflow_tracking import MLFlowTrackingService

from ..llm_service import LLMService  # existing LLM integration
from .ml_utils import plot_confusion_matrix, plot_roc_curve, plot_learning_curve, plot_feature_importance, plot_precision_recall_curve, print_classification_metrics


class MLAIClient:
    """
    ML client with rule-based and LLM-based AI suggestions.
    """

    def __init__(self, use_mlflow: bool = False, device: str = "cpu"):
        self.deep_learning_service = DeepLearningService(device=device)
        self.time_series_service = TimeSeriesService()
        self.nlp_service = NLPService()
        self.survival_service = SurvivalService()
        self.mlflow_service = MLFlowTrackingService() if use_mlflow else None
        self.llm_service = LLMService()

    # -------------------------
    # Suggest model
    # -------------------------
    def suggest_model(
        self,
        df: pd.DataFrame,
        goal: Optional[str] = None,
        mode: str = "rule"  # "rule" or "ai"
    ) -> Dict[str, Any]:
        numeric_cols = df.select_dtypes(include="number").columns.tolist()
        categorical_cols = df.select_dtypes(include="object").columns.tolist()
        text_cols = df.select_dtypes(include="object").columns.tolist()
        n_samples, n_features = df.shape

        # -------------------------
        # Rule-based suggestion
        # -------------------------
        if mode == "rule":
            if goal:
                goal = goal.lower()
                if "classification" in goal:
                    return {"service": "ml_service", "model_type": "RandomForest", "params": {}}
                elif "regression" in goal:
                    return {"service": "ml_service", "model_type": "LinearRegression", "params": {}}
                elif "nlp" in goal:
                    return {"service": "nlp_service", "model_type": "Transformer", "params": {}}
                elif "time_series" in goal:
                    return {"service": "time_series_service", "model_type": "ARIMA", "params": {}}
                elif "survival" in goal:
                    return {"service": "survival_service", "model_type": "KaplanMeier", "params": {}}
            # fallback
            return {"service": "ml_service", "model_type": "RandomForest", "params": {}}

        # -------------------------
        # LLM-based AI suggestion
        # -------------------------
        elif mode == "ai":
            prompt = f"""
            You are an expert ML engineer. 
            Dataset: {n_samples} rows, {n_features} columns
            Numeric columns: {numeric_cols}
            Categorical columns: {categorical_cols}
            Text columns: {text_cols}
            Analysis goal: {goal or 'unknown'}

            Suggest the best ML/statistical method including:
            1. service (ml_service, nlp_service, deep_learning_service, time_series_service, survival_service)
            2. model_type (e.g., RandomForest, LinearRegression, CNN, Transformer, ARIMA)
            3. optional parameters

            Respond in JSON format only.
            """
            response = self.llm_service.query(prompt)
            try:
                suggestion = eval(response)
            except Exception:
                # fallback to rule-based
                suggestion = {"service": "ml_service", "model_type": "RandomForest", "params": {}}
            return suggestion

        else:
            raise ValueError("mode must be 'rule' or 'ai'")

    # -------------------------
    # Train model
    # -------------------------
    def train_model(
        self,
        df: pd.DataFrame,
        target_col: str,
        goal: Optional[str] = None,
        mode: str = "rule",
        **kwargs
    ) -> Dict[str, Any]:
        """
        Suggests and trains a model based on the dataset.
        Automatically dispatches to the correct service.
        """
        suggestion = self.suggest_model(df, goal=goal, mode=mode)
        service_name = suggestion["service"]
        model_type = suggestion["model_type"]
        params = suggestion.get("params", {})

        X = df.drop(columns=[target_col])
        y = df[target_col]

        trained_model = None
        metrics = {}

        if service_name == "ml_service":
            from .ml_service_core import MLServiceCore  # rule-based ML methods
            ml_service = MLServiceCore()
            trained_model, metrics = ml_service.train(X, y, model_type=model_type, **params)

        elif service_name == "deep_learning_service":
            trained_model, metrics = self.deep_learning_service.train(X, y, model_type=model_type, **params)

        elif service_name == "time_series_service":
            trained_model, metrics = self.time_series_service.train(df, target_col, model_type=model_type, **params)

        elif service_name == "nlp_service":
            trained_model, metrics = self.nlp_service.train(X, y, model_type=model_type, **params)

        elif service_name == "survival_service":
            trained_model, metrics = self.survival_service.train(df, target_col, model_type=model_type, **params)

        else:
            raise ValueError(f"Unknown service: {service_name}")

        # Optional MLflow logging
        if self.mlflow_service:
            self.mlflow_service.log_model(trained_model, model_name=model_type, metrics=metrics)

        return {"trained_model": trained_model, "metrics": metrics, "suggestion": suggestion}

def train_model(
        self,
        df: pd.DataFrame,
        target_col: str,
        goal: Optional[str] = None,
        mode: str = "rule",
        **kwargs
    ) -> Dict[str, Any]:
        """
        Suggests and trains a model based on the dataset.
        Automatically dispatches to the correct service.
        """
        suggestion = self.suggest_model(df, goal=goal, mode=mode)
        service_name = suggestion["service"]
        model_type = suggestion["model_type"]
        params = suggestion.get("params", {})

        X = df.drop(columns=[target_col])
        y = df[target_col]

        trained_model = None
        metrics = {}

        if service_name == "ml_service":
            from .ml_service_core import MLServiceCore  # rule-based ML methods
            ml_service = MLServiceCore()
            trained_model, metrics = ml_service.train(X, y, model_type=model_type, **params)

        elif service_name == "deep_learning_service":
            trained_model, metrics = self.deep_learning_service.train(X, y, model_type=model_type, **params)

        elif service_name == "time_series_service":
            trained_model, metrics = self.time_series_service.train(df, target_col, model_type=model_type, **params)

        elif service_name == "nlp_service":
            trained_model, metrics = self.nlp_service.train(X, y, model_type=model_type, **params)

        elif service_name == "survival_service":
            trained_model, metrics = self.survival_service.train(df, target_col, model_type=model_type, **params)

        else:
            raise ValueError(f"Unknown service: {service_name}")

        # Optional MLflow logging
        if self.mlflow_service:
            self.mlflow_service.log_model(trained_model, model_name=model_type, metrics=metrics)

        # Plot evaluation metrics
        y_pred = trained_model.predict(X)
        if hasattr(trained_model, 'predict_proba'):  # for probability-based metrics like ROC
            y_score = trained_model.predict_proba(X)[:, 1]
        else:
            y_score = y_pred

        plot_confusion_matrix(y, y_pred)
        plot_roc_curve(y, y_score)
        plot_precision_recall_curve(y, y_score)
        plot_learning_curve(trained_model, X, y)
        plot_feature_importance(trained_model, X.columns)

        print_classification_metrics(y, y_pred)

        return {"trained_model": trained_model, "metrics": metrics, "suggestion": suggestion}