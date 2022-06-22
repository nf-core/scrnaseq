#!/usr/bin/env python
import scanpy as sc, anndata as ad, pandas as pd
from pathlib import Path
import argparse

def read_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    df.set_index("sample")

    # samplesheet may contain replicates, when it has,
    # group information from replicates and collapse with commas
    # only keep unique values using set()
    df = df.groupby(["sample"]).agg(lambda column: ",".join(set(column)))

    # return
    return(df)

# combine and write, just a simple concatenation.
def concat_h5ad(adatas):
    combined = ad.concat(adatas, label="sample", merge="unique")
    
    # return
    return(combined)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenates h5ad files and merge metadata from samplesheet")

    parser.add_argument("-i", "--input", dest="input", help="Path to samplesheet.csv")
    parser.add_argument("-o", "--out",   dest="out",   help="Output path.")

    args = vars(parser.parse_args())

    # Open samplesheet as dataframe
    df_samplesheet = read_samplesheet(args["input"])

    # find all h5ad and append to dict
    dict_of_h5ad = {
            str(path).replace("_matrix.h5ad", ""): sc.read_h5ad(path)
            for path in Path(".").rglob('*.h5ad')
    }

    # concat h5ad files
    adata = concat_h5ad(dict_of_h5ad)

    # merge with data.frame, on sample information
    adata.obs = adata.obs.join(df_samplesheet, on="sample")
    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))