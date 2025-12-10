import json
import urllib.parse
from .track_registry import TrackRegistry
from .session_builder import IGVSessionBuilder


class IGVService:

    def __init__(self):
        self.registry = TrackRegistry()
        self.builder = IGVSessionBuilder()

    def register_track(self, path, track_type, name=None, genome="hg38", metadata=None):
        return self.registry.register(
            path=path,
            track_type=track_type,
            genome=genome,
            name=name,
            metadata=metadata,
        )

    def create_session(self, genome, tracks, locus=None):
        return self.builder.build(genome, tracks, locus)

    def launch_session(self, session, mode="web"):
        if mode == "web":
            session_str = json.dumps(session)
            encoded = urllib.parse.quote(session_str)
            return f"https://igv.org/app/?session={encoded}"

        if mode == "desktop":
            path = "/tmp/igv_session.json"
            with open(path, "w") as f:
                json.dump(session, f, indent=2)
            return path

        raise ValueError("Unsupported mode")
