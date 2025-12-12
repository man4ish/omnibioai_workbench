# permissions.py
def check_permission(user, service_name):
    role_permissions = {
        "admin": ["upload", "rag_query", "report", "network_viz"],
        "researcher": ["upload", "rag_query", "report"],
        "viewer": ["rag_query"]
    }
    return service_name in role_permissions.get(user.role, [])
