import React, { useState, useEffect } from 'react';
import axios from 'axios';

const WorkflowLogs = ({ workflowId }) => {
  const [logs, setLogs] = useState([]);

  useEffect(() => {
    // Fetch logs for the selected workflow
    axios.get(`/api/workflows/${workflowId}/logs/`)
      .then(response => {
        setLogs(response.data.logs);
      })
      .catch(error => {
        console.error('There was an error fetching logs!', error);
      });
  }, [workflowId]);

  return (
    <div>
      <h2>Logs</h2>
      <ul>
        {logs.map((log, index) => (
          <li key={index}>{log}</li>
        ))}
      </ul>
    </div>
  );
};

export default WorkflowLogs;
