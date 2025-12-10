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
