"""
Module: track_registry
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides IGVTrack dataclass and TrackRegistry for managing and registering tracks
    for IGV (Integrative Genomics Viewer) sessions.

Usage:
    from omnibioai.services.track_registry import TrackRegistry

    registry = TrackRegistry()

    track = registry.register(
        path="data/sample.bam",
        track_type="alignment",
        genome="hg38",
        name="Sample BAM",
        metadata={"source": "project_X"}
    )
    print(track)
    # Output: IGVTrack(name='Sample BAM', path='data/sample.bam', type='alignment', genome='hg38', format='bam', metadata={'source': 'project_X'})

Classes:
    - IGVTrack (dataclass):
        Represents an IGV track with attributes for name, path, type, genome, format, and optional metadata.

        Attributes:
            name (str): Track name.
            path (str): File path to the track.
            type (str): Track type (e.g., alignment, variant).
            genome (str): Genome assembly (e.g., hg38).
            format (str): File format derived from the file extension.
            metadata (dict | None): Optional metadata dictionary.

    - TrackRegistry:
        Service to register IGV tracks.

        Methods:
            * register(path, track_type, genome, name=None, metadata=None) -> IGVTrack:
                Registers a track and returns an IGVTrack object.
                The track name defaults to the filename if not provided.
                The format is inferred from the file extension.
"""


from dataclasses import dataclass
from pathlib import Path

@dataclass
class IGVTrack:
    name: str
    path: str
    type: str
    genome: str
    format: str
    metadata: dict | None = None


class TrackRegistry:

    def register(self, path, track_type, genome, name=None, metadata=None):
        return IGVTrack(
            name=name or Path(path).name,
            path=path,
            type=track_type,
            genome=genome,
            format=Path(path).suffix.replace(".", ""),
            metadata=metadata,
        )
