FROM python:3.13-slim

WORKDIR /app

RUN apt update
RUN apt install cron build-essential wget libpq-dev curl -y

RUN pip install --upgrade pip

RUN pip install watchdog

COPY ./src/cellbro-db /app/cellbro-db
RUN pip install -e /app/cellbro-db

COPY ./src/cellbro-worker /app/cellbro-worker
RUN pip install -e /app/cellbro-worker