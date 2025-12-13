# tests/core_services/demo_ml_utils_save.py

from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from omnibioai.services.ml_service import ml_utils as utils
import numpy as np
import os
from datetime import datetime
# Add this at the top of your demo file
import matplotlib.pyplot as plt

# -------------------------
# Utility to save plots
# -------------------------
def save_plot(fig, name: str):
    os.makedirs("data/reports", exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = f"data/reports/{name}_{timestamp}.png"
    fig.savefig(path, bbox_inches='tight')
    print(f"[Saved] {path}")
    plt.close(fig)  # Close figure to avoid overlapping

# Patch plotting functions to return figure objects
def plot_confusion_matrix_save(*args, **kwargs):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    utils.plot_confusion_matrix(*args, **kwargs)
    save_plot(fig, "confusion_matrix")

def plot_roc_curve_save(*args, **kwargs):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    utils.plot_roc_curve(*args, **kwargs)
    save_plot(fig, "roc_curve")

def plot_precision_recall_curve_save(*args, **kwargs):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    utils.plot_precision_recall_curve(*args, **kwargs)
    save_plot(fig, "precision_recall_curve")

def plot_feature_importance_save(*args, **kwargs):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    utils.plot_feature_importance(*args, **kwargs)
    save_plot(fig, "feature_importance")

def plot_learning_curve_save(*args, **kwargs):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    utils.plot_learning_curve(*args, **kwargs)
    save_plot(fig, "learning_curve")

def plot_validation_curve_save(*args, **kwargs):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    utils.plot_validation_curve(*args, **kwargs)
    save_plot(fig, "validation_curve")


# -------------------------
# Main demo
# -------------------------
def main():
    # Load dataset
    data = load_breast_cancer()
    X, y = data.data, data.target
    feature_names = data.feature_names

    # Split train/test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

    # Train model
    model = RandomForestClassifier(n_estimators=50, random_state=42)
    model.fit(X_train, y_train)

    # Predictions
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    # -------------------------
    # Plots & metrics
    # -------------------------
    print("=== Confusion Matrix ===")
    plot_confusion_matrix_save(y_test, y_pred, labels=[0, 1])

    print("=== ROC Curve ===")
    plot_roc_curve_save(y_test, y_proba)

    print("=== Precision-Recall Curve ===")
    plot_precision_recall_curve_save(y_test, y_proba)

    print("=== Feature Importance ===")
    plot_feature_importance_save(model, feature_names, top_n=10)

    print("=== Classification Metrics ===")
    utils.print_classification_metrics(y_test, y_pred)

    print("=== Learning Curve ===")
    plot_learning_curve_save(model, X_train, y_train, cv=5, scoring='accuracy')

    print("=== Validation Curve: n_estimators ===")
    param_range = np.arange(10, 110, 20)
    plot_validation_curve_save(model, X_train, y_train, param_name='n_estimators', param_range=param_range, title='Validation Curve: n_estimators')


if __name__ == "__main__":
    main()
