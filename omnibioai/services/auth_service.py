"""
Module: auth_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides authentication utilities for generating and verifying JSON Web Tokens (JWTs) 
    for OmniBio users. Tokens encode user ID, role, and expiration time using HS256 algorithm.

Usage:
    from auth_service import AuthService
    from users.models import OmniBioUser

    # Generate token
    user = OmniBioUser.objects.get(username="john_doe")
    token = AuthService.generate_token(user, expires_minutes=120)

    # Verify token
    verified_user = AuthService.verify_token(token)
    if verified_user:
        print(f"Authenticated user: {verified_user.username}")
    else:
        print("Invalid or expired token.")

Classes:
    - AuthService:
        Static methods for handling JWT-based authentication.
        Methods:
            * generate_token(user, expires_minutes=60) -> str:
                Generates a JWT for the given user with an optional expiration time.
            * verify_token(token) -> OmniBioUser | None:
                Verifies a JWT and returns the corresponding user object if valid; otherwise, None.

Dependencies:
    - PyJWT: For encoding and decoding JSON Web Tokens.
    - Django settings: For SECRET_KEY.
    - users.models.OmniBioUser: User model for retrieving authenticated users.
"""

import jwt
from datetime import datetime, timedelta
from django.conf import settings
from users.models import OmniBioUser

SECRET_KEY = settings.SECRET_KEY
ALGORITHM = "HS256"

class AuthService:

    @staticmethod
    def generate_token(user, expires_minutes=60):
        payload = {
            "user_id": user.id,
            "role": user.role,
            "exp": datetime.utcnow() + timedelta(minutes=expires_minutes)
        }
        token = jwt.encode(payload, SECRET_KEY, algorithm=ALGORITHM)
        return token

    @staticmethod
    def verify_token(token):
        try:
            payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
            user = OmniBioUser.objects.get(id=payload["user_id"])
            return user
        except Exception as e:
            return None
