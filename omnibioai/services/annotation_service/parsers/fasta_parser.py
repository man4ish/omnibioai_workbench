"""
Module: fasta_parser
Author: Manish Kumar
Version: 1.0
Date: 2025-12-13

Description:
    Parses FASTA files to retrieve sequences by identifier or coordinates.

Classes:
    - FASTAParser:
        Methods:
            * get_sequence(identifier: str) -> str
"""

from Bio import SeqIO

class FASTAParser:
    def __init__(self, filepath):
        self.filepath = filepath
        self.sequences = {rec.id: rec.seq for rec in SeqIO.parse(filepath, "fasta")}

    def get_sequence(self, identifier: str):
        return str(self.sequences.get(identifier, ""))

