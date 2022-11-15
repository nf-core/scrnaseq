#!/usr/bin/env python
import scanpy as sc
import pandas as pd
import argparse
import os
from scipy import io
from anndata import AnnData

def _10x_h5_to_adata(mtx_h5: str, sample: str):
    adata = sc.read_10x_h5(mtx_h5)
    adata.var["gene_symbols"] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs["sample"] = sample

    # reorder columns for 10x mtx files
    adata.var = adata.var[['gene_symbols', 'feature_types', 'genome']]

    return adata


def _mtx_to_adata(
    mtx_file: str,
    barcode_file: str,
    feature_file: str,
    sample: str,
    aligner: str,
):
    adata = sc.read_mtx(mtx_file)
    if (
        aligner == "star"
    ):  # for some reason star matrix comes transposed and doesn't fit when values are appended directly
        adata = adata.transpose()

    adata.obs_names = pd.read_csv(barcode_file, header=None, sep="\t")[0].values
    adata.var_names = pd.read_csv(feature_file, header=None, sep="\t")[0].values
    adata.obs["sample"] = sample

    return adata

def input_to_adata(
    input_data: str,
    barcode_file: str,
    feature_file: str,
    sample: str,
    aligner: str,
    verbose: bool = True,
):

    if verbose:
        print("Reading in {}".format(input_data))

    if aligner == 'cellranger':
        return _10x_h5_to_adata(input_data, sample)
    else:
        return _mtx_to_adata(input_data, barcode_file, feature_file, sample, aligner)


def write_counts(
    adata: AnnData,
    txp2gene: str,
    star_index: str,
    aligner: str,
    out: str,
    verbose: bool = False,
):

    if verbose:
        print("Reading in {}".format(txp2gene))

    t2g = None
    if txp2gene:
        t2g = pd.read_table(txp2gene, header=None, names=["gene_id", "gene_symbol"], usecols=[1, 2])
    elif star_index:
        t2g = pd.read_table(f"{star_index}/geneInfo.tab", header=None, skiprows=1, names=["gene_id", "gene_symbol"], usecols=[0, 1])

    if txp2gene or star_index:
        t2g = t2g.drop_duplicates(subset='gene_id').set_index("gene_id")
        adata.var["gene_symbol"] = t2g["gene_symbol"]

    # when assigning a pandas column to a dataframe, pandas already takes care of matching the index.
    # therefore, setting the index is enough.
    pd.DataFrame(adata.obs.index).to_csv(os.path.join(out, "barcodes.tsv"), sep="\t", index=False, header=None)
    pd.DataFrame(adata.var).to_csv(os.path.join(out, "features.tsv"), sep="\t", index=True, header=None)
    io.mmwrite(os.path.join(out, "matrix.mtx"), adata.X.T, field="integer")

    if verbose:
        print("Wrote features.tsv, barcodes.tsv, and matrix.mtx files to {}".format(args["out"]))

def dump_versions():
    import pkg_resources
    with open('versions.yml', 'w') as f:
        f.write("\n".join([f"{pkg.key}: {pkg.version}" for pkg in pkg_resources.working_set]))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-i", "--input_data", dest="input_data", help="Path to either mtx or mtx h5 file.")
    parser.add_argument("-v", "--verbose", dest="verbose", help="Toggle verbose messages", default=False)
    parser.add_argument("-f", "--feature", dest="feature", help="Path to feature file.", nargs="?", const="")
    parser.add_argument("-b", "--barcode", dest="barcode", help="Path to barcode file.", nargs="?", const="")
    parser.add_argument("-s", "--sample", dest="sample", help="Sample name")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")
    parser.add_argument("-a", "--aligner", dest="aligner", help="Which aligner has been used?")
    parser.add_argument("--txp2gene", dest="txp2gene", help="Transcript to gene (t2g) file.", nargs="?", const="")
    parser.add_argument(
        "--star_index", dest="star_index", help="Star index folder containing geneInfo.tab.", nargs="?", const=""
    )

    args = vars(parser.parse_args())

    # create the directory with the sample name
    os.makedirs(os.path.dirname(args["out"]), exist_ok=True)

    print(args)
    adata = input_to_adata(
        args["input_data"], args["barcode"], args["feature"], args["sample"], args["aligner"], verbose=args["verbose"]
    )

    write_counts(adata, args["txp2gene"], args["star_index"], args["aligner"], args["sample"], verbose=args["verbose"])

    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))

    dump_versions()
