import json
from channels.generic.websocket import AsyncWebsocketConsumer

class WorkflowConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        await self.channel_layer.group_add("workflow_updates", self.channel_name)
        await self.accept()

    async def disconnect(self, close_code):
        await self.channel_layer.group_discard("workflow_updates", self.channel_name)

    # Receive updates from server
    async def workflow_update(self, event):
        await self.send(text_data=json.dumps(event["data"]))
