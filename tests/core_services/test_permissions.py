from omnibioai.services.permissions import check_permission
from types import SimpleNamespace
from omnibioai.services.logger_service import logger

# ----------------------------
# 1. Define test users
# ----------------------------
admin_user = SimpleNamespace(role="admin")
researcher_user = SimpleNamespace(role="researcher")
viewer_user = SimpleNamespace(role="viewer")
unknown_user = SimpleNamespace(role="guest")

# ----------------------------
# 2. Define services to test
# ----------------------------
services = ["upload", "rag_query", "report", "network_viz"]

# ----------------------------
# 3. Run tests
# ----------------------------
test_cases = [
    (admin_user, {"upload": True, "rag_query": True, "report": True, "network_viz": True}),
    (researcher_user, {"upload": True, "rag_query": True, "report": True, "network_viz": False}),
    (viewer_user, {"upload": False, "rag_query": True, "report": False, "network_viz": False}),
    (unknown_user, {"upload": False, "rag_query": False, "report": False, "network_viz": False}),
]

for user, expected in test_cases:
    for service_name, expected_result in expected.items():
        result = check_permission(user, service_name)
        if result == expected_result:
            logger.info(f"PASS: {user.role} -> {service_name} = {result}")
        else:
            logger.error(f"FAIL: {user.role} -> {service_name} = {result} (expected {expected_result})")

