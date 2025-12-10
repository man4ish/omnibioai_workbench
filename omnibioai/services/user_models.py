from django.contrib.auth.models import AbstractUser
from django.db import models

class OmniBioUser(AbstractUser):
    ROLE_CHOICES = (
        ('admin', 'Admin'),
        ('researcher', 'Researcher'),
        ('viewer', 'Viewer')
    )
    role = models.CharField(max_length=20, choices=ROLE_CHOICES, default='viewer')
