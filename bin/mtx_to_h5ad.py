#!/usr/bin/env python
import scanpy as sc
import pandas as pd
import argparse
import os
import scipy
from anndata import AnnData

def mtx_to_adata(
    mtx_file: str,
    barcode_file: str,
    feature_file: str,
    sample: str,
    aligner: str,
    verbose: bool = True,
):

    if verbose:
        print("Reading in {}".format(mtx_file))

    adata = sc.read_mtx(mtx_file)
    if (
        aligner == "star"
    ):  # for some reason star matrix comes transposed and doesn't fit when values are appended directly
        adata = adata.transpose()
    adata.obs_names = pd.read_csv(barcode_file, header=None, sep="\t")[0].values
    adata.var_names = pd.read_csv(feature_file, header=None, sep="\t")[0].values
    adata.obs["sample"] = sample


    return adata

def write_counts(
    adata: AnnData,
    txp2gene: str,
    out: str,
    verbose: bool = True,):

    if verbose:
        print("Reading in {}".format(txp2gene))

    try:
        os.makedirs(out)
    except FileExistsError:
        # directory already exists
        pass
    # read txp2gene file to add gene names to features file
    t2g = pd.read_table(txp2gene, header=None)
    id2name = {e[1]: e[2] for _, e in t2g.iterrows()}

    print('t2g file contents')
    print(t2g)
    # export features file with gene names
    features = pd.DataFrame()
    features["id"] = adata.var.index
    features["name"] = adata.var.index.map(id2name)
    features.to_csv(os.path.join(out, "features.tsv"), sep="\t", index=False, header=None)
    pd.DataFrame(adata.obs.index).to_csv(os.path.join(out, "barcodes.tsv"), sep="\t", index=False, header=None)
    scipy.io.mmwrite(os.path.join(out, "matrix.mtx"), adata.X.T, field="integer")

    if verbose:
        print("Wrote features.tsv, barcodes.tsv, and matrix.mtx files to {}".format(args["out"]))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-m", "--mtx", dest="mtx", help="Path to mtx file.")
    parser.add_argument("-v", "--verbose", dest="verbose", help="Toggle verbose messages", default=True)
    parser.add_argument("-f", "--feature", dest="feature", help="Path to feature file.")
    parser.add_argument("-b", "--barcode", dest="barcode", help="Path to barcode file.")
    parser.add_argument("-s", "--sample", dest="sample", help="Sample name")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")
    parser.add_argument("-a", "--aligner", dest="aligner", help="Which aligner has been used?")
    parser.add_argument("--txp2gene", dest="txp2gene", help="Transcript to gene (t2g) file.")

    args = vars(parser.parse_args())

    adata = mtx_to_adata(
        args["mtx"],
        args["barcode"],
        args["feature"],
        args["sample"],
        args["aligner"],
        verbose=args["verbose"]
    )

    write_counts(
        adata,
        args["txp2gene"],
        args["sample"],
        verbose=args["verbose"]
    )

    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))
