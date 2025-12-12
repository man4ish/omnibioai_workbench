import os
import django
import time

# ----------------------------
# Setup Django environment
# ----------------------------
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "omnibioai.settings")
django.setup()

from django.contrib.auth import get_user_model
from omnibioai.services.auth_service import AuthService
from omnibioai.services.logger_service import logger

# Use the proper user model
OmniBioUser = get_user_model()

# ----------------------------
# 1. Create test user
# ----------------------------
user, created = OmniBioUser.objects.get_or_create(
    username="testuser",
    defaults={
        "role": "tester",
        "email": "testuser@example.com",
        "password": "password123"  # only for testing
    }
)
logger.info(f"Test user created: {user.username} (created={created})")

# ----------------------------
# 2. Generate token
# ----------------------------
token = AuthService.generate_token(user, expires_minutes=5)
logger.info(f"Generated JWT token: {token}")

# ----------------------------
# 3. Verify token
# ----------------------------
verified_user = AuthService.verify_token(token)
if verified_user:
    logger.info(f"Token verified successfully for user: {verified_user.username}")
else:
    logger.error("Token verification failed")

# ----------------------------
# 4. Test expired token
# ----------------------------
expired_token = AuthService.generate_token(user, expires_minutes=0)  # expires immediately
time.sleep(1)
expired_check = AuthService.verify_token(expired_token)
if not expired_check:
    logger.info("Expired token correctly rejected")
else:
    logger.error("Expired token should not be valid")
