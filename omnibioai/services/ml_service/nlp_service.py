# ml_service/nlp_service.py
from typing import Any, Dict
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.pipeline import Pipeline

class NLPService:
    """
    NLP tasks: Naive Bayes, TF-IDF, Hugging Face
    """

    def train_naive_bayes(self, texts, labels) -> Dict[str, Any]:
        clf = Pipeline([
            ("tfidf", TfidfVectorizer()),
            ("nb", MultinomialNB())
        ])
        clf.fit(texts, labels)
        return {"model": clf, "status": "trained"}

    def train_transformer(self, texts, labels, model_name: str = "bert-base-uncased") -> Dict[str, Any]:
        # Placeholder for Hugging Face fine-tuning
        return {"status": f"Transformer {model_name} training placeholder"}
