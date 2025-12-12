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
