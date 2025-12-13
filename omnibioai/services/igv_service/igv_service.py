"""
Module: igv_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the IGVService class for managing IGV (Integrative Genomics Viewer) sessions.
    Supports registering tracks, creating IGV sessions, and launching sessions
    in either web or desktop mode.

Usage:
    from omnibioai.services.igv_service import IGVService

    igv = IGVService()

    # Register a track
    track_info = igv.register_track(
        path="data/tracks/sample.bam",
        track_type="alignment",
        name="Sample BAM",
        genome="hg38"
    )

    # Create a session
    session = igv.create_session(
        genome="hg38",
        tracks=[track_info],
        locus="chr1:100000-200000"
    )

    # Launch in web mode
    url = igv.launch_session(session, mode="web")
    print(url)

    # Launch in desktop mode
    desktop_path = igv.launch_session(session, mode="desktop")
    print(desktop_path)

Classes:
    - IGVService:
        Service for registering tracks, building IGV sessions, and launching them.

        Methods:
            * __init__():
                Initializes the IGVService with a TrackRegistry and IGVSessionBuilder.
            * register_track(path, track_type, name=None, genome="hg38", metadata=None):
                Registers a track for IGV visualization and returns track info.
            * create_session(genome, tracks, locus=None):
                Builds an IGV session with the specified genome, tracks, and optional locus.
            * launch_session(session, mode="web"):
                Launches the IGV session either in web mode (returns URL) or desktop mode
                (writes JSON file and returns file path). Raises ValueError for unsupported modes.

Dependencies:
    - json: For serializing session data.
    - urllib.parse: For URL encoding session data for web IGV.
    - TrackRegistry: Handles track registration.
    - IGVSessionBuilder: Builds IGV session objects.
"""


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
