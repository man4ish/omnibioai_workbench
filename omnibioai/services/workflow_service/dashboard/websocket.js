import React, { useState, useEffect } from 'react';

const WebSocketComponent = ({ workflowId }) => {
  const [progress, setProgress] = useState(0);

  useEffect(() => {
    const socket = new WebSocket(`ws://localhost:8000/ws/workflow/progress/${workflowId}/`);

    socket.onmessage = (event) => {
      const data = JSON.parse(event.data);
      setProgress(data.message.progress);
    };

    return () => socket.close();
  }, [workflowId]);

  return (
    <div>
      <h2>Real-time Progress</h2>
      <progress value={progress} max="100"></progress>
      <p>{progress}%</p>
    </div>
  );
};

export default WebSocketComponent;
