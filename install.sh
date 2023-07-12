#!/bin/bash

git clone https://github.com/lutrarutra/cellbro.git
git clone https://github.com/lutrarutra/scout.git
pip install -r cellbro/requirements.txt
pip install -e ./cellbro
pip install -e ./scout
cd cellbro

# Optional: Example
mkdir data && wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
python3 preprocess_example.py