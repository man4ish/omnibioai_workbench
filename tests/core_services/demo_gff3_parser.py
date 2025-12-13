from omnibioai.services.annotation_service.parsers.gff3_parser import GFF3Parser

def main():
    parser = GFF3Parser("data/sample_annotations.gff3")

    # Query a region on chr1
    results = parser.query_feature("chr1", 120, 180)
    print("Query results for chr1:120-180")
    for r in results:
        print(r)

    # Query for genes only
    gene_results = parser.query_feature("chr1", 0, 300, feature_type="gene")
    print("\nQuery for genes on chr1:0-300")
    for r in gene_results:
        print(r)

    # Query on chr2 for exons
    exon_results = parser.query_feature("chr2", 350, 360, feature_type="exon")
    print("\nQuery for exons on chr2:350-360")
    for r in exon_results:
        print(r)

if __name__ == "__main__":
    main()

