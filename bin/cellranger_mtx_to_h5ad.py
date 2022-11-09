#!/usr/bin/env python
import scanpy as sc
import argparse
import os
import pandas as pd
import scipy
from anndata import AnnData

def mtx_to_adata(mtx_h5: str, sample: str, verbose: bool = False):

    if verbose:
        print("Reading in {}".format(mtx_h5))

    adata = sc.read_10x_h5(mtx_h5)
    adata.var["gene_symbols"] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs["sample"] = sample

    return adata

def write_counts(
    adata: AnnData,
    out: str,
    verbose: bool = True,):

    pd.DataFrame(adata.var.index).to_csv(os.path.join(out, "features.tsv"), sep="\t", index=False, header=None)
    pd.DataFrame(adata.obs.index).to_csv(os.path.join(out, "barcodes.tsv"), sep="\t", index=False, header=None)
    scipy.io.mmwrite(os.path.join(out, "matrix.mtx"), adata.X.T, field="integer")

    if verbose:
        print("Wrote features.tsv, barcodes.tsv, and matrix.mtx files to {}".format(args["out"]))



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-m", "--mtx", dest="mtx", help="Path to mtx h5 file.")
    parser.add_argument("-v", "--verbose", dest="verbose", help="Toggle verbose messages", default=False)
    parser.add_argument("-s", "--sample", dest="sample", help="Sample name")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")

    args = vars(parser.parse_args())

    # create the directory with the sample name
    try:
        os.makedirs(os.path.dirname(args["out"]))
    except FileExistsError:
        # directory already exists
        pass

    adata = mtx_to_adata(args["mtx"], args["sample"], verbose=args["verbose"])

    write_counts(
        adata,
        args["sample"],
        verbose=args["verbose"]
    )

    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))
