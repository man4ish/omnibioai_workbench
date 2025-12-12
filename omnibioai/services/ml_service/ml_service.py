# manishair/services/ml_service/ml_service.py

import os
from typing import Optional, List, Dict, Any, Tuple
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    accuracy_score, roc_auc_score, f1_score, confusion_matrix,
    mean_squared_error, r2_score
)
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN

# Supervised models
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import xgboost as xgb

# Optional results directory
DEFAULT_RESULTS_DIR = os.environ.get("OMNIBIOAI_RESULTS_DIR", os.path.join("data", "ml_results"))
os.makedirs(DEFAULT_RESULTS_DIR, exist_ok=True)


class MLService:
    """
    Core ML & Statistics service.
    Provides functions for supervised & unsupervised learning, dimensionality reduction,
    and evaluation metrics.
    """

    def __init__(self, results_dir: Optional[str] = None):
        self.results_dir = results_dir or DEFAULT_RESULTS_DIR

    # -------------------------
    # Data utilities
    # -------------------------
    def _prepare_features(
        self,
        df: pd.DataFrame,
        target_col: Optional[str] = None,
        features: Optional[List[str]] = None,
        scale: bool = True
    ) -> Tuple[np.ndarray, Optional[np.ndarray], List[str]]:
        """
        Prepare feature matrix X and target y.
        Optionally scale numeric features.
        """
        if features is None:
            features = df.columns.drop(target_col) if target_col else df.columns.tolist()

        X = df[features].copy()
        numeric_cols = X.select_dtypes(include=np.number).columns.tolist()
        categorical_cols = X.select_dtypes(include="object").columns.tolist()

        # Encode categorical columns
        for col in categorical_cols:
            le = LabelEncoder()
            X[col] = le.fit_transform(X[col].astype(str))

        # Scale numeric columns
        if scale and numeric_cols:
            scaler = StandardScaler()
            X[numeric_cols] = scaler.fit_transform(X[numeric_cols])

        y = df[target_col].values if target_col else None
        return X.values, y, features

    # -------------------------
    # Supervised Learning
    # -------------------------
    def train_classification(
        self,
        df: pd.DataFrame,
        target_col: str,
        features: Optional[List[str]] = None,
        model_type: str = "random_forest",
        test_size: float = 0.2,
        random_state: int = 42,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Train a classification model and return metrics and predictions.
        Supported model_types: 'logistic_regression', 'random_forest', 'xgboost'
        """
        X, y, used_features = self._prepare_features(df, target_col, features)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)

        # Select model
        if model_type == "logistic_regression":
            model = LogisticRegression(**kwargs)
        elif model_type == "random_forest":
            model = RandomForestClassifier(**kwargs)
        elif model_type == "xgboost":
            model = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', **kwargs)
        else:
            raise ValueError(f"Unsupported model_type: {model_type}")

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)[:, 1] if hasattr(model, "predict_proba") else None

        metrics = {
            "accuracy": accuracy_score(y_test, y_pred),
            "f1_score": f1_score(y_test, y_pred, average="weighted"),
            "roc_auc": roc_auc_score(y_test, y_proba) if y_proba is not None else None,
            "confusion_matrix": confusion_matrix(y_test, y_pred).tolist()
        }

        return {
            "model": model,
            "metrics": metrics,
            "predictions": y_pred,
            "probabilities": y_proba,
            "features": used_features
        }

    def train_regression(
        self,
        df: pd.DataFrame,
        target_col: str,
        features: Optional[List[str]] = None,
        model_type: str = "linear_regression",
        test_size: float = 0.2,
        random_state: int = 42,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Train a regression model and return metrics and predictions.
        Supported model_types: 'linear_regression', 'random_forest', 'xgboost'
        """
        X, y, used_features = self._prepare_features(df, target_col, features)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)

        # Select model
        if model_type == "linear_regression":
            model = LinearRegression(**kwargs)
        elif model_type == "random_forest":
            model = RandomForestRegressor(**kwargs)
        elif model_type == "xgboost":
            model = xgb.XGBRegressor(**kwargs)
        else:
            raise ValueError(f"Unsupported model_type: {model_type}")

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        metrics = {
            "mse": mean_squared_error(y_test, y_pred),
            "rmse": np.sqrt(mean_squared_error(y_test, y_pred)),
            "r2": r2_score(y_test, y_pred)
        }

        return {
            "model": model,
            "metrics": metrics,
            "predictions": y_pred,
            "features": used_features
        }

    # -------------------------
    # Unsupervised Learning
    # -------------------------
    def train_clustering(
        self,
        df: pd.DataFrame,
        features: Optional[List[str]] = None,
        method: str = "kmeans",
        n_clusters: int = 3,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Train an unsupervised clustering model.
        Supported methods: 'kmeans', 'dbscan'
        """
        X, _, used_features = self._prepare_features(df, features=features)

        if method == "kmeans":
            model = KMeans(n_clusters=n_clusters, **kwargs)
            labels = model.fit_predict(X)
        elif method == "dbscan":
            model = DBSCAN(**kwargs)
            labels = model.fit_predict(X)
        else:
            raise ValueError(f"Unsupported clustering method: {method}")

        return {
            "model": model,
            "labels": labels,
            "features": used_features
        }

    # -------------------------
    # Dimensionality Reduction
    # -------------------------
    def apply_pca(
        self,
        df: pd.DataFrame,
        features: Optional[List[str]] = None,
        n_components: int = 2,
        scale: bool = True
    ) -> Dict[str, Any]:
        """
        Apply PCA to reduce dimensionality.
        Returns transformed components and explained variance.
        """
        X, _, used_features = self._prepare_features(df, features=features, scale=scale)
        pca = PCA(n_components=n_components)
        components = pca.fit_transform(X)

        return {
            "pca_model": pca,
            "components": components,
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "features": used_features
        }
