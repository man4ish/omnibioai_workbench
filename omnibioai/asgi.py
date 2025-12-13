from channels.routing import ProtocolTypeRouter, URLRouter
from channels.auth import AuthMiddlewareStack
from django.urls import path
from omnibioai.services.workflow_service.dashboard.websocket import WorkflowProgressConsumer

application = ProtocolTypeRouter({
    "http": get_asgi_application(),
    "websocket": AuthMiddlewareStack(
        URLRouter([
            path("ws/workflow/progress/<int:workflow_id>/", WorkflowProgressConsumer.as_asgi()),
        ])
    ),
})
