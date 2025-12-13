import subprocess
import os

class SnakemakeAdapter:
    """
    Adapter to run Snakemake workflows with progress reporting
    """

    def __init__(self, snakemake_bin="snakemake"):
        self.snakemake_bin = snakemake_bin

    def run(self, entrypoint, work_dir=None, config=None, progress=False, params=None):
        """
        Run Snakemake workflow. Yields (progress, log) if progress=True.

        Parameters
        ----------
        entrypoint : str
            Path to Snakefile
        work_dir : str
            Working directory
        config : dict
            Config dictionary for Snakemake
        progress : bool
            Yield step progress
        params : dict
            Extra params like {"cores": 1}
        """
        cmd = [self.snakemake_bin, "--snakefile", entrypoint, "--printshellcmds", "--rerun-incomplete"]

        # Set cores
        cores = 1
        if params and "cores" in params:
            cores = params["cores"]
        cmd.extend(["--cores", str(cores)])

        # Config
        if config:
            for k, v in config.items():
                cmd.extend(["--config", f"{k}={v}"])

        if work_dir:
            os.makedirs(work_dir, exist_ok=True)
            cmd.extend(["--directory", work_dir])

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        if progress:
            step = 0
            for line in process.stdout:
                line = line.strip()
                step += 1
                yield min(step, 100), line
        else:
            output, _ = process.communicate()
            if process.returncode != 0:
                raise RuntimeError(f"Snakemake workflow failed: {output}")
            yield 100, output


