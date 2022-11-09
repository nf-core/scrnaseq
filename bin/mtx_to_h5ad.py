#!/usr/bin/env python
import scanpy as sc
import pandas as pd
import argparse


def mtx_to_adata(
    mtx_file: str,
    barcode_file: str,
    feature_file: str,
    sample: str,
    aligner: str,
    verbose: bool = False,
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-m", "--mtx", dest="mtx", help="Path to mtx file.")
    parser.add_argument("-v", "--verbose", dest="verbose", help="Toggle verbose messages", default=False)
    parser.add_argument("-f", "--feature", dest="feature", help="Path to feature file.")
    parser.add_argument("-b", "--barcode", dest="barcode", help="Path to barcode file.")
    parser.add_argument("-s", "--sample", dest="sample", help="Sample name")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")
    parser.add_argument("-a", "--aligner", dest="aligner", help="Which aligner has been used?")

    args = vars(parser.parse_args())

    adata = mtx_to_adata(
        args["mtx"],
        args["barcode"],
        args["feature"],
        args["sample"],
        args["aligner"],
        verbose=args["verbose"],
    )

    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))
