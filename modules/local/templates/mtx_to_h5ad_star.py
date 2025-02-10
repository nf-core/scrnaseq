#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
import argparse
import anndata
from anndata import AnnData, concat
import platform
from scipy.sparse import csr_matrix

def _mtx_to_adata(
    input: str,
    sample: str,
):
    adata = sc.read_10x_mtx(input)
    adata.obs["sample"] = sample
    adata.layers["count"] = adata.X.copy()

    velocyto_dir = f"velocyto_{input}"
    if os.path.exists(velocyto_dir):
        barcodes = os.path.join(velocyto_dir, "barcodes.tsv.gz")
        features = os.path.join(velocyto_dir, "features.tsv.gz")

        for matrix in ["ambiguous", "spliced", "unspliced"]:
            adata_state = sc.read_mtx(os.path.join(velocyto_dir, f"{matrix}.mtx.gz")).T

            adata_state.obs_names = pd.read_csv(barcodes, header=None, sep="\\t")[0].values
            adata_state.var_names = pd.read_csv(features, header=None, sep="\\t")[0].values

            missing_obs = adata.obs_names[~adata.obs_names.isin(adata_state.obs_names)]
            adata_missing = AnnData(
                X=csr_matrix((len(missing_obs), adata.shape[1])),
                obs=pd.DataFrame(index=missing_obs),
                var=adata_state.var
            )
            adata_state = concat([adata_state, adata_missing], join="outer")
            adata_state = adata_state[adata.obs_names, adata.var["gene_ids"]].copy()

            adata.layers[matrix] = adata_state.X

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
    adata.var["gene_symbol"] = adata.var.index
    adata.var['gene_versions'] = adata.var["gene_ids"]
    adata.var.index = adata.var['gene_versions'].str.split('.').str[0].values
    adata.var_names_make_unique()  # in case user does not use ensembl references, names might not be unique

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
    input_data="${meta.input_type}",
    output="${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}"
)

# dump versions
dump_versions()
