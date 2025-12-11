from xml.etree.ElementTree import Element, SubElement, tostring
from pathlib import Path

class IGVSessionBuilder:
    """
    Creates IGV session XML files from BAM/VCF/FASTA tracks.
    """
    
    def __init__(self, genome="hg38"):
        self.genome = genome
        self.tracks = []

    def add_track(self, path: str, track_type="alignment"):
        self.tracks.append({"path": path, "type": track_type})

    def build(self, output_path: str):
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

        session = Element("Session", genome=self.genome, version="8")

        resources = SubElement(session, "Resources")
        for t in self.tracks:
            SubElement(resources, "Resource", path=t["path"])

        xml = tostring(session, encoding="unicode")
        Path(output_path).write_text(xml)

        return xml
