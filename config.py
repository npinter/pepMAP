import os


class Config:
    SECRET_KEY = os.urandom(24)
    PEPMAP_STORE_DIR = "storage"
    PEPMAP_STORE_TTL_SECONDS = 1800  # 30 minutes
    START_SCHEDULER = True
    CLEAR_STORE_ON_START = True
    APP_HOST = "0.0.0.0"
    APP_PORT = 7007
    DEBUG = False
