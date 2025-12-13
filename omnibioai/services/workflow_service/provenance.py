import mysql.connector
import json
from datetime import datetime

class WorkflowProvenance:
    def __init__(self, host="localhost", user="root", password="", database="omnibioai"):
        self.conn = mysql.connector.connect(
            host=host,
            user=user,
            password=password,
            database=database
        )
        self.cursor = self.conn.cursor(dictionary=True)

    def create_run(self, workflow_name, engine, entrypoint, parameters):
        self.start_time = datetime.now()
        self.run_id = None
        self.logs = []
        self.outputs = []

        sql = """
            INSERT INTO workflow_runs (workflow_name, engine, entrypoint, parameters, start_time, state)
            VALUES (%s, %s, %s, %s, %s, %s)
        """
        self.cursor.execute(sql, (
            workflow_name,
            engine,
            entrypoint,
            json.dumps(parameters),
            self.start_time,
            "RUNNING"
        ))
        self.conn.commit()
        self.run_id = self.cursor.lastrowid
        return self.run_id

    def update_progress(self, progress, log_line):
        self.logs.append(log_line)
        sql = "UPDATE workflow_runs SET logs=%s WHERE id=%s"
        self.cursor.execute(sql, (json.dumps(self.logs), self.run_id))
        self.conn.commit()

    def complete_run(self, outputs):
        self.outputs = outputs
        end_time = datetime.now()
        sql = """
            UPDATE workflow_runs
            SET state=%s, end_time=%s, outputs=%s, logs=%s
            WHERE id=%s
        """
        self.cursor.execute(sql, (
            "COMPLETED",
            end_time,
            json.dumps(self.outputs),
            json.dumps(self.logs),
            self.run_id
        ))
        self.conn.commit()

    def fail_run(self, error_message):
        self.logs.append(error_message)
        end_time = datetime.now()
        sql = """
            UPDATE workflow_runs
            SET state=%s, end_time=%s, logs=%s
            WHERE id=%s
        """
        self.cursor.execute(sql, (
            "FAILED",
            end_time,
            json.dumps(self.logs),
            self.run_id
        ))
        self.conn.commit()

