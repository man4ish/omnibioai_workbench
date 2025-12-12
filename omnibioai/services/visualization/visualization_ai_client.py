# omnibioai/services/visualization_ai_client.py

from typing import Optional, Dict, Any, Tuple
import pandas as pd
import json

from .plot_dispatch import plot_dispatch
from .llm_service import LLMService  # assumes you have a service for LLM calls


class VisualizationAIClient:
    """
    AI-assisted visualization client.
    Suggests plots based on dataset characteristics or user goal,
    and creates the plot using `plot_dispatch`.
    """

    def __init__(self, results_dir: Optional[str] = None, llm: Optional[LLMService] = None):
        self.results_dir = results_dir
        self.llm = llm or LLMService()  # default LLM service

    # -------------------------
    # 1️⃣ Suggest plot using LLM
    # -------------------------
    def suggest_plot(
        self,
        df: pd.DataFrame,
        goal: Optional[str] = None
    ) -> Tuple[str, Dict[str, Any]]:
        """
        Suggest a plot type and parameters using AI.
        """
        # summarize dataset
        numeric_cols = df.select_dtypes(include="number").columns.tolist()
        categorical_cols = df.select_dtypes(include="object").columns.tolist()
        n_rows, n_cols = df.shape
        df_summary = {
            "n_rows": n_rows,
            "n_cols": n_cols,
            "numeric_cols": numeric_cols,
            "categorical_cols": categorical_cols,
            "goal": goal
        }

        # prompt the AI
        prompt = f"""
        You are a data visualization expert.
        I have a dataset with the following properties:
        {json.dumps(df_summary)}

        Suggest the best plot type for this data and any parameters needed
        (e.g., x_col, y_col, sample_column, group_col). 
        Return a JSON object like:
        {{
            "plot_type": "<plot_type>",
            "plot_params": {{}}
        }}
        """
        # call LLM
        response_text = self.llm.ask(prompt)
        try:
            response_json = json.loads(response_text)
            plot_type = response_json.get("plot_type")
            plot_params = response_json.get("plot_params", {})
        except Exception:
            # fallback to rules
            plot_type, plot_params = self._fallback_suggest(df, goal)

        return plot_type, plot_params

    # -------------------------
    # fallback rules (old heuristic)
    # -------------------------
    def _fallback_suggest(self, df: pd.DataFrame, goal: Optional[str] = None) -> Tuple[str, Dict[str, Any]]:
        numeric_cols = df.select_dtypes(include="number").columns.tolist()
        categorical_cols = df.select_dtypes(include="object").columns.tolist()
        plot_params = {}

        if goal == "variance" or goal == "dimensionality":
            if len(numeric_cols) >= 2:
                plot_type = "pca"
                plot_params["sample_column"] = categorical_cols[0] if categorical_cols else None
                return plot_type, plot_params

        if goal == "correlation":
            if len(numeric_cols) >= 2:
                plot_type = "heatmap"
                plot_params["row_names"] = df.index.name
                plot_params["col_names"] = None
                return "heatmap", plot_params

        if len(numeric_cols) >= 2 and len(categorical_cols) >= 1:
            plot_type = "boxplot"
            plot_params["value_col"] = numeric_cols[0]
            plot_params["group_col"] = categorical_cols[0]
        elif len(numeric_cols) >= 2:
            plot_type = "scatter"
            plot_params["x_col"] = numeric_cols[0]
            plot_params["y_col"] = numeric_cols[1]
        elif len(numeric_cols) == 1 and len(categorical_cols) >= 1:
            plot_type = "bar"
            plot_params["x_col"] = categorical_cols[0]
            plot_params["y_col"] = numeric_cols[0]
        else:
            plot_type = "histogram"
            plot_params["value_col"] = numeric_cols[0] if numeric_cols else None

        return plot_type, plot_params

    # -------------------------
    # 2️⃣ Generate plot
    # -------------------------
    def create_plot(
        self,
        df: pd.DataFrame,
        plot_type: Optional[str] = None,
        save_path_prefix: Optional[str] = None,
        return_types: Optional[list] = None,
        goal: Optional[str] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Create a plot using either AI-suggested or specified plot_type.
        """
        if plot_type is None:
            plot_type, suggested_params = self.suggest_plot(df, goal=goal)
            suggested_params.update(kwargs)
        else:
            suggested_params = kwargs

        return plot_dispatch(
            plot_type=plot_type,
            df=df,
            save_path_prefix=save_path_prefix,
            return_types=return_types,
            results_dir=self.results_dir,
            **suggested_params
        )
