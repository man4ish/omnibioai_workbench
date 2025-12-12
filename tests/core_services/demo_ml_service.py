# tests/core_services/demo_all_services.py

import pandas as pd
import numpy as np
from sklearn.datasets import load_breast_cancer, load_diabetes
from ml_service.ml_ai_client import MLAIClient

def demo_classification(ml_client):
    print("\n=== Classification Demo ===")
    data = load_breast_cancer(as_frame=True)
    df = data.frame
    target_col = 'target'

    # Rule-based
    print("\n--- Rule-based ---")
    suggestion_rule = ml_client.suggest_model(df, goal="classification", mode="rule")
    print("Suggestion:", suggestion_rule)
    result_rule = ml_client.train_model(df, target_col, goal="classification", mode="rule")
    print("Metrics:", result_rule["metrics"])

    # AI-based
    print("\n--- AI-based ---")
    suggestion_ai = ml_client.suggest_model(df, goal="classification", mode="ai")
    print("Suggestion:", suggestion_ai)
    result_ai = ml_client.train_model(df, target_col, goal="classification", mode="ai")
    print("Metrics:", result_ai["metrics"])


def demo_regression(ml_client):
    print("\n=== Regression Demo ===")
    data = load_diabetes(as_frame=True)
    df = data.frame
    target_col = 'target'

    # Rule-based
    print("\n--- Rule-based ---")
    suggestion_rule = ml_client.suggest_model(df, goal="regression", mode="rule")
    print("Suggestion:", suggestion_rule)
    result_rule = ml_client.train_model(df, target_col, goal="regression", mode="rule")
    print("Metrics:", result_rule["metrics"])

    # AI-based
    print("\n--- AI-based ---")
    suggestion_ai = ml_client.suggest_model(df, goal="regression", mode="ai")
    print("Suggestion:", suggestion_ai)
    result_ai = ml_client.train_model(df, target_col, goal="regression", mode="ai")
    print("Metrics:", result_ai["metrics"])


def demo_deep_learning(ml_client):
    print("\n=== Deep Learning Demo ===")
    X = np.random.rand(10, 1, 28, 28)  # Example: 10 samples, 1 channel, 28x28 images
    y = np.random.randint(0, 10, size=(10,))

    print("\n--- CNN ---")
    res = ml_client.deep_learning_service.train_cnn(X, y)
    print(res)

    print("\n--- LSTM ---")
    res = ml_client.deep_learning_service.train_lstm(X, y)
    print(res)

    print("\n--- Transformer ---")
    res = ml_client.deep_learning_service.train_transformer(["sample text"], [0])
    print(res)


def demo_nlp(ml_client):
    print("\n=== NLP Demo ===")
    texts = ["I love ML", "AI is amazing", "I dislike bugs"]
    labels = [1, 1, 0]

    print("\n--- Naive Bayes ---")
    res = ml_client.nlp_service.train_naive_bayes(texts, labels)
    print(res)

    print("\n--- Transformer ---")
    res = ml_client.nlp_service.train_transformer(texts, labels)
    print(res)


def demo_time_series(ml_client):
    print("\n=== Time Series Demo ===")
    series = pd.Series(np.random.randn(20).cumsum())

    print("\n--- ARIMA ---")
    res = ml_client.time_series_service.train_arima(series)
    print(res)

    print("\n--- SARIMA ---")
    res = ml_client.time_series_service.train_sarima(series)
    print(res)

    print("\n--- Bayesian ---")
    res = ml_client.time_series_service.train_bayesian(series)
    print(res)


def demo_survival(ml_client):
    print("\n=== Survival Analysis Demo ===")
    durations = pd.Series([5,6,6,2,4])
    events = pd.Series([1,0,0,1,1])
    df = pd.DataFrame({
        "duration": durations,
        "event": events,
        "age": [50, 60, 65, 45, 70]
    })

    print("\n--- Kaplan-Meier ---")
    res = ml_client.survival_service.fit_kaplan_meier(durations, events)
    print(res)

    print("\n--- Cox Regression ---")
    res = ml_client.survival_service.fit_cox_regression(df, duration_col="duration", event_col="event")
    print(res)


if __name__ == "__main__":
    np.random.seed(42)
    ml_client = MLAIClient(use_mlflow=False, device="cpu")

    demo_classification(ml_client)
    demo_regression(ml_client)
    demo_deep_learning(ml_client)
    demo_nlp(ml_client)
    demo_time_series(ml_client)
    demo_survival(ml_client)

    print("\n=== All Services Demo Completed ===")
