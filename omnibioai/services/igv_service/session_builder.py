
"""
Module: session_builder
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the IGVSessionBuilder class to construct IGV (Integrative Genomics Viewer)
    session objects from genome information and registered tracks. Supports optional
    locus specification for focused visualization.

Usage:
    from omnibioai.services.session_builder import IGVSessionBuilder

    builder = IGVSessionBuilder()

    # Example track object list
    tracks = [
        type("Track", (), {"name": "Sample BAM", "path": "data/sample.bam", "type": "alignment", "format": "bam"})()
    ]

    session = builder.build(genome="hg38", tracks=tracks, locus="chr1:100000-200000")
    print(session)
    # Output: {'genome': 'hg38', 'tracks': [{'name': 'Sample BAM', 'url': 'data/sample.bam', 'type': 'alignment', 'format': 'bam'}], 'locus': 'chr1:100000-200000'}

Classes:
    - IGVSessionBuilder:
        Builds IGV session dictionaries for visualization in web or desktop IGV.

        Methods:
            * build(genome, tracks, locus=None) -> dict:
                Constructs a session object including genome, tracks, and optional locus.
                Each track should have attributes: name, path, type, format.
"""

class IGVSessionBuilder:

    def build(self, genome, tracks, locus=None):
        session = {
            "genome": genome,
            "tracks": [
                {
                    "name": t.name,
                    "url": t.path,
                    "type": t.type,
                    "format": t.format,
                }
                for t in tracks
            ],
        }

        if locus:
            session["locus"] = locus

        return session
