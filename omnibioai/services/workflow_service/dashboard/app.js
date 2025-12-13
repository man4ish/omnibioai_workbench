import React, { useState, useEffect } from 'react';
import axios from 'axios';
import WorkflowProgress from './components/WorkflowProgress';
import WorkflowLogs from './components/WorkflowLogs';

function App() {
  const [workflows, setWorkflows] = useState([]);
  const [selectedWorkflow, setSelectedWorkflow] = useState(null);

  useEffect(() => {
    // Fetch the list of workflows from the backend
    axios.get('/api/workflows/')
      .then(response => {
        setWorkflows(response.data);
      })
      .catch(error => {
        console.error('There was an error fetching workflows!', error);
      });
  }, []);

  const handleSelectWorkflow = (workflow) => {
    setSelectedWorkflow(workflow);
  };

  return (
    <div className="App">
      <h1>OmnibioAI Workflow Dashboard</h1>
      
      <div>
        <h2>Workflows</h2>
        <ul>
          {workflows.map(workflow => (
            <li key={workflow.id}>
              <button onClick={() => handleSelectWorkflow(workflow)}>
                {workflow.workflow_name} ({workflow.state})
              </button>
            </li>
          ))}
        </ul>
      </div>

      {selectedWorkflow && (
        <div>
          <WorkflowProgress workflowId={selectedWorkflow.id} />
          <WorkflowLogs workflowId={selectedWorkflow.id} />
        </div>
      )}
    </div>
  );
}

export default App;
