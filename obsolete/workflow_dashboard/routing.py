from django.urls import path
from .consumers import WorkflowConsumer

websocket_urlpatterns = [
    path("ws/workflows/", WorkflowConsumer.as_asgi()),
]
