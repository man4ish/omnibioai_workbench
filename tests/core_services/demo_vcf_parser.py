# File: tests/core_services/demo_vcf_parser.py

"""
Demo for VCFParser
"""

from omnibioai.services.annotation_service.parsers.vcf_parser import VCFParser

def main():
    # Path to your tabix-indexed VCF file
    vcf_file = "data/reference/fake.vcf.gz"

    parser = VCFParser(vcf_file)

    # Single variant query
    chrom = "1"
    pos = 1234567
    ref = "A"
    alt = "T"

    result = parser.query_variant(chrom, pos, ref, alt)
    print("Query result for variant 1:1234567 A>T")
    print(result)

    # Another variant
    chrom2 = "2"
    pos2 = 7654321
    ref2 = "G"
    alt2 = "C"

    result2 = parser.query_variant(chrom2, pos2, ref2, alt2)
    print("\nQuery result for variant 2:7654321 G>C")
    print(result2)

if __name__ == "__main__":
    main()
