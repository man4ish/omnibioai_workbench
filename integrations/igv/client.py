import socket
from .exceptions import IGVConnectionError, IGVCommandError

class IGVClient:
    """
    Controls IGV via its socket API.
    """

    def __init__(self, host="localhost", port=60151, timeout=5):
        self.host = host
        self.port = port
        self.timeout = timeout

    # -------------------------------
    # Connect to IGV
    # -------------------------------
    def _connect(self):
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.settimeout(self.timeout)
            s.connect((self.host, self.port))
            return s
        except Exception as e:
            raise IGVConnectionError(
                f"Failed to connect to IGV at {self.host}:{self.port} - {str(e)}"
            )

    # -------------------------------
    # Send IGV Command
    # -------------------------------
    def send(self, command: str) -> str:
        sock = self._connect()
        try:
            sock.sendall((command + "\n").encode())
            response = sock.recv(4096).decode().strip()
            sock.close()

            if response.startswith("ERROR"):
                raise IGVCommandError(f"IGV returned error: {response}")

            return response

        except Exception as e:
            raise IGVCommandError(f"IGV command failed: {command} -> {str(e)}")

    # -------------------------------
    # Convenience wrappers
    # -------------------------------
    def load(self, path: str):
        return self.send(f"load {path}")

    def set_genome(self, genome: str):
        return self.send(f"genome {genome}")

    def goto(self, locus: str):
        return self.send(f"goto {locus}")

    def snapshot(self, filename: str):
        return self.send(f"snapshot {filename}")

    def snapshot_directory(self, path: str):
        return self.send(f"snapshotDirectory {path}")
