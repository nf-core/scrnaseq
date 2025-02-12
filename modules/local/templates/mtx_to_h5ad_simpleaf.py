#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
import argparse
import anndata
from anndata import AnnData
import platform
import json

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
            "pandas": pd.__version__,
            "anndata": anndata.__version__,
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
    simpleaf_h5ad_path = f"{input_data}/alevin/quants.h5ad"

    # the simpleaf quant module exports an h5ad file.
    adata = sc.read_h5ad(simpleaf_h5ad_path)
    adata.obs_names = adata.obs['barcodes'].values
    adata.obs["sample"] = sample

    # standard format
    # index are gene IDs and symbols are a column
    adata.var['gene_versions'] = adata.var['gene_id']
    adata.var.index = adata.var['gene_versions'].str.split('.').str[0].values
    adata.var_names_make_unique()

    # sort adata column- and row- wise to avoid positional differences
    adata = adata[adata.obs_names.sort_values(), adata.var_names.sort_values()]

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
    input_data="${inputs}",
    output="${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}"
)

# dump versions
dump_versions()
