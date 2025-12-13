"""
Demo: TSVParser
File: tests/core_services/demo_tsv_parser.py
Author: Manish Kumar
Date: 2025-12-13
"""

import os
from omnibioai.services.annotation_service.parsers.tsv_parser import TSVParser

def create_sample_tsv():
    os.makedirs("data", exist_ok=True)
    sample_data = """gene\tchrom\tstart\tend\tscore
geneA\tchr1\t100\t200\t0.5
geneB\tchr1\t150\t250\t0.8
geneC\tchr2\t300\t400\t0.7
geneD\tchr2\t350\t450\t0.9
"""
    with open("data/sample_annotations.tsv", "w") as f:
        f.write(sample_data)
    print("Sample TSV file created at data/sample_annotations.tsv")

def main():
    # Create sample TSV
    create_sample_tsv()

    # Initialize parser
    parser = TSVParser("data/sample_annotations.tsv")

    # Query by column filter
    print("Query results for chrom=chr1:")
    results = parser.query({"chrom": "chr1"})
    for r in results:
        print(r)

    print("\nQuery results for score=0.7:")
    results = parser.query({"score": 0.7})
    for r in results:
        print(r)

if __name__ == "__main__":
    main()

