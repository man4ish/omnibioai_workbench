def format_region(chrom: str, pos: int, window: int = 100):
    return f"{chrom}:{pos-window}-{pos+window}"
