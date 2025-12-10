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
