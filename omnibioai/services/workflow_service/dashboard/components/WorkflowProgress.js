import React, { useState, useEffect } from 'react';
import axios from 'axios';

const WorkflowProgress = ({ workflowId }) => {
  const [progress, setProgress] = useState(0);

  useEffect(() => {
    const interval = setInterval(() => {
      axios.get(`/api/workflows/${workflowId}/`)
        .then(response => {
          const workflow = response.data;
          setProgress(workflow.state === 'RUNNING' ? workflow.progress : 100);
        })
        .catch(error => {
          console.error('There was an error fetching progress!', error);
        });
    }, 5000); // Update every 5 seconds

    return () => clearInterval(interval);
  }, [workflowId]);

  return (
    <div>
      <h2>Progress</h2>
      <progress value={progress} max="100"></progress>
      <p>{progress}%</p>
    </div>
  );
};

export default WorkflowProgress;
