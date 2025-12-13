# demo_nlp_service.py
import sys
import os

# Ensure root directory is on Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from omnibioai.services.ml_service.nlp_service import NLPService
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

def save_confusion_matrix(y_true, y_pred, filename="nb_confusion_matrix.png"):
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Naive Bayes Confusion Matrix")
    save_path = os.path.join("data", "reports", f"{filename}")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path)
    plt.close()
    print(f"[Saved] {save_path}")

def main():
    texts = ["I love this!", "I hate that!"]
    labels = ["positive", "negative"]

    nlp_service = NLPService()
    result = nlp_service.train_naive_bayes(texts, labels)
    clf = result["model"]

    # Predictions
    y_pred = clf.predict(texts)
    
    # Classification report
    print("=== Naive Bayes Classification Report ===")
    print(classification_report(labels, y_pred))

    # Save confusion matrix
    save_confusion_matrix(labels, y_pred)

    # Transformer placeholder
    transformer_result = nlp_service.train_transformer(texts, labels)
    print("\n=== Transformer Placeholder ===")
    print(transformer_result["status"])

if __name__ == "__main__":
    main()
