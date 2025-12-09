import pandas as pd
from sklearn.linear_model import LogisticRegression

# Example: Dummy model
def make_prediction(df):
    df = df.select_dtypes(include='number')  # Dummy cleanup
    model = LogisticRegression()
    X = df.iloc[:, :-1]
    y = df.iloc[:, -1]
    model.fit(X, y)
    preds = model.predict(X)
    df['Prediction'] = preds
    return df
