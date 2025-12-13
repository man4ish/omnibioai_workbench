import subprocess
import os

class WDLAdapter:
    """
    Adapter to run WDL workflows using Cromwell with progress reporting
    """

    def __init__(self, cromwell_jar_path="cromwell.jar", java_bin="java"):
        self.cromwell_jar_path = cromwell_jar_path
        self.java_bin = java_bin

    def run_workflow(self, wdl_file, inputs_json=None, options_json=None, work_dir=None, progress=False):
        """
        Run a WDL workflow with optional inputs/options JSON.
        Yields (progress, log) if progress=True
        """
        cmd = [self.java_bin, "-jar", self.cromwell_jar_path, "run", wdl_file]

        if inputs_json:
            cmd.extend(["-i", inputs_json])
        if options_json:
            cmd.extend(["-o", options_json])

        if work_dir:
            os.makedirs(work_dir, exist_ok=True)
            os.chdir(work_dir)

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
                raise RuntimeError(f"WDL workflow failed: {output}")
            yield 100, output
