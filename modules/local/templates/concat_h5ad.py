#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc, anndata as ad, pandas as pd
from pathlib import Path


def read_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    df.set_index("sample")

    # samplesheet may contain replicates, when it has,
    # group information from replicates and collapse with commas
    # only keep unique values using set()
    df = df.groupby(["sample"]).agg(lambda column: ",".join(set(column.astype(str)))

    return df


if __name__ == "__main__":

    # Open samplesheet as dataframe
    df_samplesheet = read_samplesheet("${samplesheet}")

    # find all h5ad and append to dict
    dict_of_h5ad = {str(path).replace("_matrix.h5ad", ""): sc.read_h5ad(path) for path in Path(".").rglob("*.h5ad")}

    # concat h5ad files
    adata = ad.concat(dict_of_h5ad, label="sample", merge="unique", index_unique="_")

    # merge with data.frame, on sample information
    adata.obs = adata.obs.join(df_samplesheet, on="sample", how="left").astype(str)
    adata.write_h5ad("combined_${meta.input_type}_matrix.h5ad")

    print("Wrote h5ad file to {}".format("combined_${meta.input_type}_matrix.h5ad"))
