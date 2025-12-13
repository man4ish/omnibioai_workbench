"""
Module: permissions
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides a simple role-based access control (RBAC) utility.
    Determines whether a user has permission to access a given service
    based on their role.

Usage:
    from omnibioai.services.permissions import check_permission

    if check_permission(user, "report"):
        print("Access granted")
    else:
        print("Access denied")

Functions:
    - check_permission(user, service_name) -> bool:
        Checks whether the given user has access to the specified service.
        Returns True if allowed, False otherwise.

Roles and Permissions:
    - admin: ["upload", "rag_query", "report", "network_viz"]
    - researcher: ["upload", "rag_query", "report"]
    - viewer: ["rag_query"]

Dependencies:
    None (expects `user` to have a `role` attribute).
"""

def check_permission(user, service_name):
    role_permissions = {
        "admin": ["upload", "rag_query", "report", "network_viz"],
        "researcher": ["upload", "rag_query", "report"],
        "viewer": ["rag_query"]
    }
    return service_name in role_permissions.get(user.role, [])
