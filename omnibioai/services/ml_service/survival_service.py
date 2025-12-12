# ml_service/survival_service.py
from typing import Any, Dict
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter

class SurvivalService:
    """
    Survival analysis: Kaplan-Meier, Cox regression
    """

    def fit_kaplan_meier(self, durations: pd.Series, event_observed: pd.Series) -> Dict[str, Any]:
        kmf = KaplanMeierFitter()
        kmf.fit(durations, event_observed)
        return {"model": kmf, "status": "fitted"}

    def fit_cox_regression(self, df: pd.DataFrame, duration_col: str, event_col: str) -> Dict[str, Any]:
        cph = CoxPHFitter()
        cph.fit(df, duration_col=duration_col, event_col=event_col)
        return {"model": cph, "status": "fitted"}
