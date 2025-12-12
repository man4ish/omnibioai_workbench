from django.contrib.auth.models import AbstractUser
from django.db import models

ROLE_CHOICES = [("viewer", "Viewer"), ("tester", "Tester")]

class OmniBioUser(AbstractUser):
    role = models.CharField(max_length=20, choices=ROLE_CHOICES, default="viewer")
