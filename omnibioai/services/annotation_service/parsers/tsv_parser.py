"""
Module: tsv_parser
Author: Manish Kumar
Version: 1.0
Date: 2025-12-13

Description:
    Generic parser for tab-delimited or CSV annotation files.
    Supports filtering by any column values.

Classes:
    - TSVParser:
        Methods:
            * query(filter_dict: dict) -> list:
                Returns rows matching filter criteria.
"""

import pandas as pd

class TSVParser:
    def __init__(self, filepath, sep='\t'):
        self.filepath = filepath
        self.df = pd.read_csv(filepath, sep=sep)

    def query(self, filter_dict: dict):
        df = self.df.copy()
        for col, val in filter_dict.items():
            df = df[df[col] == val]
        return df.to_dict(orient="records")

