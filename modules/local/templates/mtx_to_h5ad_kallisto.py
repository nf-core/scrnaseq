#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
import argparse
from anndata import AnnData
import platform
import glob

def _mtx_to_adata(
    input: str,
    sample: str,
    t2g: str
):

    adata = sc.read_mtx(glob.glob(f"{input}/*.mtx")[0])
    adata.obs_names = pd.read_csv(glob.glob(f"{input}/*.barcodes.txt")[0], header=None, sep="\\t")[0].values
    adata.var_names = pd.read_csv(glob.glob(f"{input}/*.genes.txt")[0], header=None, sep="\\t")[0].values
    adata.obs["sample"] = sample

    txp2gene = pd.read_table(glob.glob(f"{t2g}")[0], header=None, names=["gene_id", "gene_symbol"], usecols=[1, 2])
    txp2gene = txp2gene.drop_duplicates(subset="gene_id").set_index("gene_id")
    adata.var = pd.merge(adata.var, txp2gene, left_index=True, right_index=True)

    return adata


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.
    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.
    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

def dump_versions():
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "pandas": pd.__version__
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))

def input_to_adata(
    input_data: str,
    output: str,
    sample: str,
    t2g: str
):
    print(f"Reading in {input_data}")

    # open main data
    adata = _mtx_to_adata(input_data, sample, t2g)

    # standard format
    # index are gene IDs and symbols are a column
    adata.var['gene_versions'] = adata.var.index
    adata.var['gene_ids'] = adata.var['gene_versions'].str.split('.').str[0]
    adata.var.index = adata.var["gene_ids"].values
    adata.var = adata.var.drop("gene_ids", axis=1)
    adata.var_names_make_unique()

    # write results
    adata.write_h5ad(f"{output}", compression="gzip")
    print(f"Wrote h5ad file to {output}")

    # dump versions
    dump_versions()

    return adata

#
# Run main script
#

# create the directory with the sample name
os.makedirs("${meta.id}", exist_ok=True)

# input_type comes from NF module
adata = input_to_adata(
    input_data="${inputs}",
    output="${meta.id}/${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}",
    t2g="${txp2gene}"
)
