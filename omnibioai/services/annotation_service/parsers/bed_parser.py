"""
Module: bed_parser
Author: Manish Kumar
Version: 1.0
Date: 2025-12-13

Description:
    Parses BED files for genomic region annotation.
    Provides query capabilities for overlapping regions.

Classes:
    - BEDParser:
        Methods:
            * query_region(chrom: str, start: int, end: int) -> list:
                Returns features overlapping the given region.
"""

import pandas as pd

class BEDParser:
    def __init__(self, filepath):
        self.filepath = filepath
        self.df = pd.read_csv(filepath, sep='\t', header=None,
                              names=["chrom", "start", "end", "name", "score", "strand"])

    def query_region(self, chrom: str, start: int, end: int):
        overlapping = self.df[(self.df["chrom"] == chrom) &
                              (self.df["start"] <= end) &
                              (self.df["end"] >= start)]
        return overlapping.to_dict(orient="records")

