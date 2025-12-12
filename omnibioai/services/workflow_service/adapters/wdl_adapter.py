import subprocess
import os

class WDLAdapter:
    """
    Adapter to run WDL workflows using Cromwell
    """

    def __init__(self, cromwell_jar_path="cromwell.jar", java_bin="java"):
        self.cromwell_jar_path = cromwell_jar_path
        self.java_bin = java_bin

    def run_workflow(self, wdl_file, inputs_json=None, options_json=None, work_dir=None):
        """
        Run a WDL workflow with optional inputs and options JSON files
        """
        cmd = [self.java_bin, "-jar", self.cromwell_jar_path, "run", wdl_file]

        if inputs_json:
            cmd.extend(["-i", inputs_json])
        if options_json:
            cmd.extend(["-o", options_json])

        if work_dir:
            os.makedirs(work_dir, exist_ok=True)
            # Cromwell output dir can be set using options_json or CWD
            os.chdir(work_dir)

        process = subprocess.run(cmd, capture_output=True, text=True)
        if process.returncode != 0:
            raise RuntimeError(f"WDL workflow failed: {process.stderr}")

        return process.stdout
