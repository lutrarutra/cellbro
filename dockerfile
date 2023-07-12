FROM ubuntu:latest

RUN apt update && apt install -y \
    python3 \
    python3-pip \
    git \
    wget

RUN pip3 install --upgrade pip
RUN git clone https://github.com/lutrarutra/cellbro.git
COPY data/data.py cellbro/data.py
RUN git clone git@github.com:lutrarutra/scout.git
RUN pip3 install -r cellbro/requirements.txt
RUN pip3 install -e ./cellbro
RUN pip3 install -e ./scout
RUN mkdir cellbro/data && wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O cellbro/data/pbmc3k_filtered_gene_bc_matrices.tar.gz
RUN python3 
RUN tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/

CMD ["python3", "cellbro/run.py" data/adata.py]
