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
    matrix: str,
    barcodes: str,
    features: str,
    t2g: str,
    sample: str
):

    adata = sc.read_mtx(matrix)
    adata.obs_names = pd.read_csv(barcodes, header=None, sep="\\t")[0].values
    adata.var_names = pd.read_csv(features, header=None, sep="\\t")[0].values
    adata.obs["sample"] = sample

    txp2gene = pd.read_table(glob.glob(f"{t2g}")[0], header=None, names=["gene_id", "gene_symbol"], usecols=[1, 2])
    txp2gene = txp2gene.drop_duplicates(subset="gene_id").set_index("gene_id")
    adata.var = adata.var.join(txp2gene, how="left")

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
    matrix: str,
    barcodes: str,
    features: str,
    t2g: str,
    output: str,
    sample: str,
):
    print(f"Reading in {matrix}")

    # open main data
    adata = _mtx_to_adata(matrix=matrix, barcodes=barcodes, features=features, sample=sample, t2g=t2g)

    # standard format
    # index are gene IDs and symbols are a column
    adata.var['gene_versions'] = adata.var.index
    adata.var.index = adata.var['gene_versions'].str.split('.').str[0]
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
if "${params.kb_workflow}" == "standard":
    input_to_adata(
        matrix=glob.glob("${inputs}/*.mtx")[0],
        barcodes=glob.glob("${inputs}/*.barcodes.txt")[0],
        features=glob.glob("${inputs}/*.genes.txt")[0],
        output="${meta.id}/${meta.id}_${meta.input_type}_matrix.h5ad",
        sample="${meta.id}",
        t2g="${txp2gene}"
    )

else:
    for type in ['spliced', 'unspliced']:
        input_to_adata(
            matrix=glob.glob("${inputs}/" + f"{type}*.mtx")[0],
            barcodes=glob.glob("${inputs}/" + f"{type}*.barcodes.txt")[0],
            features=glob.glob("${inputs}/" + f"{type}*.genes.txt")[0],
            output="${meta.id}/${meta.id}_${meta.input_type}" + f"_{type}_matrix.h5ad",
            sample="${meta.id}",
            t2g="${txp2gene}"
        )

# dump versions
dump_versions()
