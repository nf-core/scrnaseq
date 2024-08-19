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
    mtx_file: str,
    barcode_file: str,
    feature_file: str,
    sample: str,
):
    adata = sc.read_mtx(mtx_file)
    adata = adata.transpose()
    adata.obs_names = pd.read_csv(barcode_file, header=None, sep="\t")[0].values
    adata.var_names = pd.read_csv(feature_file, header=None, sep="\t")[0].values
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
    barcode_file: str,
    feature_file: str,
    sample: str,
    txp2gene: str,
):
    print(f"Reading in {input_data}")

    # open main data
    adata = _mtx_to_adata(input_data, barcode_file, feature_file, sample)

    # open gene information
    print(f"Reading in {txp2gene}")
    t2g = pd.read_table(f"{txp2gene}", header=None, skiprows=1, names=["gene_id", "gene_symbol"], usecols=[0, 1])
    t2g = t2g.drop_duplicates(subset="gene_id").set_index("gene_id")

    # standard format
    adata.var.index = t2g.index.values.tolist()
    adata.var["gene_symbol"] = t2g["gene_symbol"]

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
    input_data="${input_type}/matrix.mtx.gz",
    output="${meta.id}/${meta.id}_${input_type}_matrix.h5ad",
    barcode_file="${input_type}/barcodes.tsv.gz",
    feature_file="${input_type}/features.tsv.gz",
    sample="${meta.id}",
    txp2gene="${star_index}/geneInfo.tab"
)
