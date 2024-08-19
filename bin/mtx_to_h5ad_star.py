#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
import argparse
from scipy import io
from anndata import AnnData

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


def input_to_adata(
    input_data: str,
    barcode_file: str,
    feature_file: str,
    sample: str,
    txp2gene: str,
):
    print("Reading in {}".format(input_data))

    # open main data
    adata = _mtx_to_adata(input_data, barcode_file, feature_file, sample)

    # open gene information
    print("Reading in {}".format(txp2gene))
    t2g = pd.read_table("{}".format(txp2gene), header=None, skiprows=1, names=["gene_id", "gene_symbol"], usecols=[0, 1])
    t2g = t2g.drop_duplicates(subset="gene_id").set_index("gene_id")

    # standard format
    adata.var.index = t2g.index.values.tolist()
    adata.var["gene_symbol"] = t2g["gene_symbol"]

    return adata

def dump_versions(task_process):
    import pkg_resources

    with open("versions.yml", "w") as f:
        f.write(f"{task_process}:\n\t")
        f.write("\n\t".join([f"{pkg.key}: {pkg.version}" for pkg in pkg_resources.working_set]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-i", "--input_data", dest="input_data", help="Path to either mtx or mtx h5 file.")
    parser.add_argument("-f", "--feature", dest="feature", help="Path to feature file.", nargs="?", const="")
    parser.add_argument("-b", "--barcode", dest="barcode", help="Path to barcode file.", nargs="?", const="")
    parser.add_argument("-s", "--sample", dest="sample", help="Sample name")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")
    parser.add_argument("--task_process", dest="task_process", help="Task process name.")
    parser.add_argument("--txp2gene", dest="txp2gene", help="Star index folder containing geneInfo.tab.", nargs="?", const="")

    args = vars(parser.parse_args())

    # create the directory with the sample name
    os.makedirs(os.path.dirname(args["out"]), exist_ok=True)

    adata = input_to_adata(
        input_data=args["input_data"],
        barcode_file=args["barcode"],
        feature_file=args["feature"],
        sample=args["sample"],
        txp2gene=args["txp2gene"]
    )

    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))

    dump_versions(task_process=args["task_process"])
