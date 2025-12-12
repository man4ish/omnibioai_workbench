# tests/core_services/demo_annotation_service.py
import pandas as pd
from omnibioai.services.annotation_service import AnnotationService

def main():
    # Initialize the service
    annotation_service = AnnotationService()

    # -------------------------
    # Single variant annotation
    # -------------------------
    chrom, pos, ref, alt = "1", 12345, "A", "G"
    annotation = annotation_service.annotate_variant(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        dataset_name="gnomAD"
    )
    print("Single variant annotation:")
    print(annotation)

    # -------------------------
    # Batch variant annotation
    # -------------------------
    df = pd.DataFrame([
        {"chrom": "1", "pos": 12345, "ref": "A", "alt": "G"},
        {"chrom": "2", "pos": 54321, "ref": "C", "alt": "T"},
    ])

    batch_results = annotation_service.annotate_variants_batch(df, dataset_name="gnomAD")
    print("\nBatch variant annotation results:")
    print(batch_results)

if __name__ == "__main__":
    main()
