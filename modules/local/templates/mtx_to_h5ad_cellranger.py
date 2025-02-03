#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os
import scanpy as sc
import pandas as pd
import argparse
import anndata
from anndata import AnnData
import platform
import glob

os.environ["NUMBA_CACHE_DIR"] = "."

def _mtx_to_adata(
    input: str,
    sample: str,
):

    adata = sc.read_10x_h5(input)
    adata.var["gene_symbols"] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs["sample"] = sample

    # reorder columns for 10x mtx files
    adata.var = adata.var[["gene_symbols", "feature_types", "genome"]]

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
    adata = _mtx_to_adata(input_data, sample)

    # standard format
    # index are gene IDs and symbols are a column
    adata.var['gene_versions'] = adata.var.index
    adata.var.index = adata.var['gene_versions'].str.split('.').str[0].values
    adata.var_names_make_unique()

    # write results
    adata.write_h5ad(f"{output}")
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
    input_data=glob.glob("*${meta.input_type}_feature_bc_matrix.h5")[0], # cellrangermulti has 'sample_' as prefix
    output="${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}"
)
