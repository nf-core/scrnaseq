#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
import argparse
from anndata import AnnData
import platform

def _mtx_to_adata(
    input: str,
    sample: str,
):

    adata = sc.read_mtx(f"{input}/quants_mat.mtx")
    adata.obs_names = pd.read_csv(f"{input}/quants_mat_rows.txt", header=None, sep="\\t")[0].values
    adata.var_names = pd.read_csv(f"{input}/quants_mat_cols.txt", header=None, sep="\\t")[0].values
    adata.obs["sample"] = sample

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
):
    print(f"Reading in {input_data}")

    # open main data
    adata = _mtx_to_adata(input_data, sample)

    # standard format
    # index are gene IDs and symbols are a column
    # TODO: how to get gene_symbols for alevin?
    adata.var['gene_versions'] = adata.var.index
    adata.var.index = adata.var['gene_versions'].str.split('.').str[0].values
    adata.var_names_make_unique()

    # write results
    adata.write_h5ad(f"{output}")
    print(f"Wrote h5ad file to {output}")

#
# Run main script
#

# create the directory with the sample name
os.makedirs("${meta.id}", exist_ok=True)

# input_type comes from NF module
input_to_adata(
    input_data="${meta.id}_alevin_results/af_quant/alevin/",
    output="${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}"
)

# dump versions
dump_versions()
