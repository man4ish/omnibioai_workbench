"""
Module: constants
Author: Manish Kumar
Description:
    Defines global constants used throughout the OmniBioAI Workbench project, including:
    - Default chromosomes for genomic analysis
    - Standard VCF file columns
    - Default number of PCA components for analysis
"""

DEFAULT_CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
VCF_COLUMNS = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]
DEFAULT_PCA_COMPONENTS = 2
