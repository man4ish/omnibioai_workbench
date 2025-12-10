import jwt
from datetime import datetime, timedelta
from django.conf import settings
from .user_models import OmniBioUser

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
