#!/usr/bin/env python
import scanpy as sc
import argparse

def mtx_to_adata( mtx_dir: str, sample: str, verbose: bool = False ):

    if verbose:
        print("Reading in {}".format(mtx_dir))

    adata = sc.read_10x_mtx(mtx_dir)
    adata.obs["sample"] = sample

    return adata


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-m", "--mtx",     dest="mtx",     help="Path to mtx directory."                 )
    parser.add_argument("-v", "--verbose", dest="verbose", help="Toggle verbose messages", default=False )
    parser.add_argument("-s", "--sample",  dest="sample",  help="Sample name"                            )
    parser.add_argument("-o", "--out",     dest="out",     help="Output path."                           )

    args = vars(parser.parse_args())

    adata = mtx_to_adata(args["mtx"], args["sample"], verbose=args["verbose"])

    adata.write_h5ad(args["out"], compression="gzip")

    print("Wrote h5ad file to {}".format(args["out"]))
