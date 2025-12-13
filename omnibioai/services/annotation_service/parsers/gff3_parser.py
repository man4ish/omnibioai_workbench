"""
Module: gff3_parser
Author: Manish Kumar
Version: 1.0
Date: 2025-12-13

Description:
    Parses GFF3 files for genomic feature annotation.
    Supports querying by region, feature type, or gene.

Classes:
    - GFF3Parser:
        Methods:
            * query_feature(chrom: str, start: int, end: int, feature_type: str = None)
"""

import pandas as pd

class GFF3Parser:
    def __init__(self, filepath):
        self.filepath = filepath
        self.df = pd.read_csv(filepath, sep='\t', comment='#', header=None,
                              names=["seqid", "source", "type", "start", "end", 
                                     "score", "strand", "phase", "attributes"])

    def query_feature(self, chrom: str, start: int, end: int, feature_type: str = None):
        df = self.df[(self.df["seqid"] == chrom) &
                     (self.df["start"] <= end) &
                     (self.df["end"] >= start)]
        if feature_type:
            df = df[df["type"] == feature_type]
        return df.to_dict(orient="records")

