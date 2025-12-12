# ml_service/ml_utils.py
from typing import Optional
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_curve,
    auc,
    confusion_matrix,
    RocCurveDisplay,
    ConfusionMatrixDisplay,
)
from sklearn.model_selection import learning_curve
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.model_selection import validation_curve
import numpy as np


def plot_confusion_matrix(y_true, y_pred, labels=None, normalize: bool = False, title: Optional[str] = None):
    cm = confusion_matrix(y_true, y_pred, labels=labels, normalize='true' if normalize else None)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=labels)
    disp.plot(cmap=plt.cm.Blues)
    plt.title(title or "Confusion Matrix")
    plt.tight_layout()
    plt.show()


def plot_roc_curve(y_true, y_score, pos_label=1, title: Optional[str] = None):
    fpr, tpr, _ = roc_curve(y_true, y_score, pos_label=pos_label)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title or 'Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.show()


def plot_learning_curve(estimator, X, y, cv=5, scoring=None, title: Optional[str] = None):
    train_sizes, train_scores, test_scores = learning_curve(estimator, X, y, cv=cv, scoring=scoring)
    train_mean = np.mean(train_scores, axis=1)
    train_std = np.std(train_scores, axis=1)
    test_mean = np.mean(test_scores, axis=1)
    test_std = np.std(test_scores, axis=1)

    plt.figure()
    plt.plot(train_sizes, train_mean, 'o-', color='r', label="Training score")
    plt.plot(train_sizes, test_mean, 'o-', color='g', label="Cross-validation score")
    plt.fill_between(train_sizes, train_mean - train_std, train_mean + train_std, alpha=0.1, color='r')
    plt.fill_between(train_sizes, test_mean - test_std, test_mean + test_std, alpha=0.1, color='g')
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    plt.title(title or "Learning Curve")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()


def plot_feature_importance(model, feature_names=None, top_n: int = 20, title: Optional[str] = None):
    if hasattr(model, "feature_importances_"):
        importances = model.feature_importances_
    elif hasattr(model, "coef_"):
        importances = np.abs(model.coef_).flatten()
    else:
        raise ValueError("Model does not have feature_importances_ or coef_ attribute")

    if feature_names is None:
        feature_names = [f"f{i}" for i in range(len(importances))]

    indices = np.argsort(importances)[::-1][:top_n]
    plt.figure(figsize=(10, 6))
    sns.barplot(x=importances[indices], y=np.array(feature_names)[indices])
    plt.title(title or "Top Feature Importances")
    plt.tight_layout()
    plt.show()



def plot_precision_recall_curve(y_true, y_score, pos_label=1, title: Optional[str] = None):
    precision, recall, _ = precision_recall_curve(y_true, y_score, pos_label=pos_label)
    average_precision = average_precision_score(y_true, y_score, pos_label=pos_label)
    
    plt.figure()
    plt.plot(recall, precision, color='b', lw=2, label=f'Precision-Recall curve (AP = {average_precision:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(title or 'Precision-Recall Curve')
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.show()


def print_classification_metrics(y_true, y_pred, average='binary'):
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, average=average)
    recall = recall_score(y_true, y_pred, average=average)
    f1 = f1_score(y_true, y_pred, average=average)
    
    print(f"Accuracy: {accuracy:.2f}")
    print(f"Precision: {precision:.2f}")
    print(f"Recall: {recall:.2f}")
    print(f"F1-Score: {f1:.2f}")



def plot_validation_curve(estimator, X, y, param_name, param_range, title: Optional[str] = None):
    train_scores, test_scores = validation_curve(estimator, X, y, param_name=param_name, param_range=param_range, cv=5)
    train_mean = np.mean(train_scores, axis=1)
    test_mean = np.mean(test_scores, axis=1)
    train_std = np.std(train_scores, axis=1)
    test_std = np.std(test_scores, axis=1)
    
    plt.figure()
    plt.plot(param_range, train_mean, label='Training score', color='r')
    plt.plot(param_range, test_mean, label='Cross-validation score', color='g')
    plt.fill_between(param_range, train_mean - train_std, train_mean + train_std, alpha=0.1, color='r')
    plt.fill_between(param_range, test_mean - test_std, test_mean + test_std, alpha=0.1, color='g')
    plt.xlabel(param_name)
    plt.ylabel("Score")
    plt.title(title or f"Validation Curve for {param_name}")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()
