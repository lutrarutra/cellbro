from pydantic_settings import BaseSettings

import pytz

class Settings(BaseSettings):
    SECRET_KEY: str = ""
    TZ: str = "UTC"
    SESSION_EXPIRE_SECONDS: int = 60 * 60 * 24 * 7  # 7 days
    JWT_ALGORITHM: str = "HS256"

    # Database
    POSTGRES_USER: str = "cellbro"
    POSTGRES_PASSWORD: str = "password"
    POSTGRES_HOST: str = "localhost"
    POSTGRES_DB: str = "cellbro_db"
    POSTGRES_PORT: int = 5432
    
    # Redis
    REDIS_HOST: str = "localhost"
    REDIS_PORT: int = 6379

    # read from environment variables
    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"

    @property
    def TIMEZONE(self) -> pytz.BaseTzInfo:
        return pytz.timezone(self.TZ)

settings = Settings()