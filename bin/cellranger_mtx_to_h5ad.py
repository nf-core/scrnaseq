#!/usr/bin/env python
import scanpy as sc
import argparse
import os
import pandas as pd
from scipy import io
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
    txp2gene: str,
    out: str,
    verbose: bool = True,):

    features = pd.DataFrame()
    features["id"] = adata.var.index

    # if txp2gene file is available enrich features file with gene names
    if txp2gene:
        t2g = pd.read_table(f"{txp2gene}/star/geneInfo.tab", header=None, skiprows=1)
        print(t2g)
        id2name = {e[0]: e[1] for _, e in t2g.iterrows()}
        print(id2name)
        features["name"] = adata.var.index.map(id2name)

    features.to_csv(os.path.join(out, "features.tsv"), sep="\t", index=False, header=None)
    pd.DataFrame(adata.obs.index).to_csv(os.path.join(out, "barcodes.tsv"), sep="\t", index=False, header=None)
    io.mmwrite(os.path.join(out, "matrix.mtx"), adata.X.T, field="integer")

    if verbose:
        print("Wrote features.tsv, barcodes.tsv, and matrix.mtx files to {}".format(args["out"]))



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-m", "--mtx", dest="mtx", help="Path to mtx h5 file.")
    parser.add_argument("-v", "--verbose", dest="verbose", help="Toggle verbose messages", default=False)
    parser.add_argument("-s", "--sample", dest="sample", help="Sample name")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")
    parser.add_argument("--txp2gene", dest="txp2gene", help="Transcript to gene (t2g) file.", nargs='?', const='')

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
        args["txp2gene"],
        args["sample"],
        verbose=args["verbose"]
    )

    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))
