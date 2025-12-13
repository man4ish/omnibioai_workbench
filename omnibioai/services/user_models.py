"""
Module: users.models
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Defines the OmniBioUser model, extending Django's AbstractUser to include a role field.
    Supports role-based access control (RBAC) for the OmniBioAI application.

Usage:
    from users.models import OmniBioUser

    # Creating a new user
    user = OmniBioUser.objects.create_user(
        username="john_doe",
        password="secure_password",
        role="tester"
    )

    # Accessing user role
    print(user.role)  # Output: "tester"

Classes:
    - OmniBioUser(AbstractUser):
        Custom user model extending Django's AbstractUser with an additional 'role' field.

        Fields:
            * role (CharField): User role, with choices "viewer" or "tester". Defaults to "viewer".

Meta:
    - app_label: "omnibioai" to specify the Django app this model belongs to.

Constants:
    - ROLE_CHOICES: List of tuples defining allowed roles for users.
"""


from django.contrib.auth.models import AbstractUser
from django.db import models

ROLE_CHOICES = [("viewer", "Viewer"), ("tester", "Tester")]

class OmniBioUser(AbstractUser):
    role = models.CharField(max_length=20, choices=ROLE_CHOICES, default="viewer")

    class Meta:
        app_label = "omnibioai"  # set this to the Django app name
