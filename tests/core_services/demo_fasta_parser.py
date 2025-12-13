"""
Demo for FASTAParser
"""

from omnibioai.services.annotation_service.parsers.fasta_parser import FASTAParser

def main():
    fasta_file = "data/sample_sequences.fasta"

    # Initialize parser
    parser = FASTAParser(fasta_file)

    # Query existing sequence
    seq_id_1 = "seq1"
    seq1 = parser.get_sequence(seq_id_1)
    print(f"Sequence for {seq_id_1}: {seq1}")

    # Query another existing sequence
    seq_id_2 = "seq2"
    seq2 = parser.get_sequence(seq_id_2)
    print(f"Sequence for {seq_id_2}: {seq2}")

    # Query a non-existing sequence
    missing_id = "seqX"
    seq_missing = parser.get_sequence(missing_id)
    if not seq_missing:
        print(f"Sequence for {missing_id} not found.")
    else:
        print(f"Sequence for {missing_id}: {seq_missing}")

if __name__ == "__main__":
    main()
