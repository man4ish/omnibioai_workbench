import os
from channels.routing import ProtocolTypeRouter, URLRouter
from channels.auth import AuthMiddlewareStack
from django.core.asgi import get_asgi_application
import plugins.workflow_dashboard.routing as dashboard_routing

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "omnibioai.settings")

application = ProtocolTypeRouter({
    "http": get_asgi_application(),
    "websocket": AuthMiddlewareStack(
        URLRouter(dashboard_routing.websocket_urlpatterns)
    ),
})
