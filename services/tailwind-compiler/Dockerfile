FROM node:25-slim
WORKDIR /app
RUN apt-get update && apt-get install -y procps curl && \
    curl -sLO https://github.com/tailwindlabs/tailwindcss/releases/latest/download/tailwindcss-linux-arm64 && \
    chmod +x tailwindcss-linux-arm64 && mv tailwindcss-linux-arm64 /usr/local/bin/tailwindcss && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
    
COPY ./services/tailwind-compiler/package*.json ./
RUN npm install
COPY ./services/tailwind-compiler/tailwind.config.js ./
RUN mkdir -p /app/static/style/compiled