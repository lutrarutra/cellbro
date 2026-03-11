from passlib.context import CryptContext
import jwt
import datetime as dt

from .config import settings

pwd_context = CryptContext(schemes=["argon2"], deprecated="auto")

def hash_password(password: str) -> str:
    return pwd_context.hash(password)

def verify_password(plain_password: str, hashed_password: str) -> bool:
    return pwd_context.verify(plain_password, hashed_password)

def create_password_reset_token(user_id: int, valid_minutes: int = 60 * 24) -> str:
    expire = dt.datetime.now() + dt.timedelta(minutes=valid_minutes)
    payload = {
        "user_id": user_id,
        "exp": expire,
        "action": "password_reset"
    }
    if not settings.SECRET_KEY:
        raise ValueError("SECRET_KEY is not set in settings.")
    return jwt.encode(payload, settings.SECRET_KEY, algorithm=settings.JWT_ALGORITHM)

def verify_password_reset_token(token: str) -> int | None:
    try:
        payload = jwt.decode(token, settings.SECRET_KEY, algorithms=[settings.JWT_ALGORITHM])
        if payload.get("action") != "password_reset":
            return None
        return payload.get("user_id")
    except jwt.ExpiredSignatureError:
        return None
    except jwt.InvalidTokenError:
        return None