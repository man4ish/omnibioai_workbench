"""
Module: vcf_parser
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the VCFParser class for querying indexed VCF files using pysam.TabixFile.
    Supports retrieval of variant annotations for given chromosome, position, reference, and alternate alleles.

Usage:
    from omnibioai.services.parsers.vcf_parser import VCFParser

    parser = VCFParser("data/reference/gnomad.vcf.gz")

    variant_info = parser.query_variant("1", 1234567, "A", "T")
    print(variant_info)
    # Output: {'chrom': '1', 'pos': 1234567, 'ref': 'A', 'alt': 'T', 'info': {'AF': '0.01', ...}}

Classes:
    - VCFParser:
        Provides variant-level queries for tabix-indexed VCF files.
        
        Methods:
            * __init__(vcf_path: str):
                Initializes the parser and opens the VCF file with pysam.TabixFile.
            * query_variant(chrom: str, pos: int, ref: str, alt: str) -> dict:
                Queries a specific variant and returns annotation as a dictionary.
                If the variant is not found, 'info' will be None.

Dependencies:
    - pysam: For querying tabix-indexed VCF files.
    - typing.Dict: For type hinting the return value.
"""


import pysam
from typing import Dict

class VCFParser:
    """
    Query indexed VCF files using pysam.TabixFile
    """
    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.tabix = pysam.TabixFile(vcf_path)

    def query_variant(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """
        Query a variant in the VCF file and return annotation as dict
        """
        try:
            records = self.tabix.fetch(chrom, pos-1, pos)  # 0-based indexing
            for record in records:
                cols = record.strip().split('\t')
                v_ref = cols[3]
                v_alt = cols[4]
                if v_ref == ref and v_alt == alt:
                    info = {k: v for k,v in (x.split('=') for x in cols[7].split(';') if '=' in x)}
                    return {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "info": info
                    }
            return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "info": None}
        except ValueError:
            return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "info": None}
