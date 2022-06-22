#!/usr/bin/env python
import scanpy as sc, anndata as ad, pandas as pd
from pathlib import Path
import argparse

def read_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    return(df)

# combine and write
# combination without inner or out join, just a simple concatenation.
def concat_h5ad(outfile):
    combined = ad.concat(list_of_h5ad)
    return(combined)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenates h5ad files and merge metadata from samplesheet")

    parser.add_argument("-i", "--input", dest="input", help="Path to samplesheet.csv")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")

    args = vars(parser.parse_args())

    # how to merge this on adata.obs?
    df_samplesheet = read_samplesheet(args["input"])

    # find all and append to list
    list_of_h5ad = [sc.read_h5ad(path) for path in Path(".").rglob('*.h5ad')]

    # concat and write
    adata = concat_h5ad(args["out"])
    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))