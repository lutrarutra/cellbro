FROM python:latest

RUN apt update && apt install --no-install-recommends -y \
    build-essential \
    gcc \
    autoconf \
    automake \
    cmake \
    flex \
    libc-dev \
    libigraph-dev \
    git \
    wget \
    curl

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN pip install --upgrade pip
RUN pip install --upgrade setuptools==58.2.0 wheel
RUN git clone https://github.com/lutrarutra/cellbro.git
RUN git clone https://github.com/lutrarutra/scout.git
RUN pip install -r cellbro/requirements.txt
RUN pip install -e ./cellbro
RUN pip install -e ./scout
RUN cd cellbro
RUN mkdir data && wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
RUN tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
RUN python3 preprocess_example.py

CMD ["python3", "cellbro/run.py", "--host 0.0.0.0", data/adata.py]
