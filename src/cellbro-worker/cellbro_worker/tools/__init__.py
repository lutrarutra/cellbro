import os
from redis import Redis

worker_redis = Redis(host="redis-cache", port=int(os.environ["REDIS_PORT"]), db=5, decode_responses=True)

from . import dataset, wrapper