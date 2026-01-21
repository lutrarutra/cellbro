FROM python:3.13-slim
WORKDIR /app
RUN apt update && apt install build-essential libpq-dev curl postgresql-client -y && rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip

COPY ./src/cellbro-db /app/cellbro-db
COPY ./src/cellbro-server /app/cellbro-server
COPY ./src/cellbro-worker /app/cellbro-worker
RUN pip install -e /app/cellbro-db
RUN pip install -e /app/cellbro-server
RUN pip install -e /app/cellbro-worker
RUN mkdir -p /static/style/compiled

CMD ["uvicorn", "cellbro_server.main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]