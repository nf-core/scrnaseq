#!/usr/bin/env python
import argparse
import sys
import os
from typing import NamedTuple
from xmlrpc.client import Boolean
import pandas as pd
import scanpy as sc
from emptydrops import find_nonambient_barcodes
from emptydrops.matrix import CountMatrix


"""
Requires 3 files within single directory
"legacy"
-barcodes.tsv
-genes.tsv
-matrix.mtx
OR
"v3"
-barcodes.tsv.gz
-features.tsv.gz
-matrix.mtx.gz
"""

### TODO: parser may change for alevin-fry
class GeneralCountMatrix(CountMatrix):
    """
    Inherits from CountMatrix, adds parsing method for alevin
    """

    @staticmethod
    def from_alevin(genome_dir: str):
        """
        Args:
            genome_dir: path to directory containing matrix, gene/feature, and barcode equivalent files from alevin
        TODO: Add parser for alevin to cell ranger type output
        """
        from vpolo.alevin import parser

        ### alevin https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py
        ### kallisto utilizes from_legacy_mtx
        ### star utilizes from_v3_mtx
        ### Cellranger is accounted for by from_legacy_mtx and from_v3_mtx


def find_candidates(matrix: str, filter_params: list = None):
    """
    Uses emptydrops python implementation to identify nonambient barcodes (empty)
    Returns namedtuple of candidate barcode indices and associated statistics
    """
    if len(filter_params) < 3:
        sys.exit("Not enough filtering parameters included, please include 3 and try again")

    return find_nonambient_barcodes(
        matrix,  # Full expression matrix
        matrix.bcs,  # (iterable of str): Strings of initially-called cell barcodes
        min_umi_frac_of_median=filter_params[0],
        min_umis_nonambient=filter_params[1],
        max_adj_pvalue=filter_params[2],
    )


def filter_candidates(data_obj: CountMatrix, candidates_list: NamedTuple, outfile: str, legacy: Boolean):
    """
    Filters matrix and barcodes based on candidate list index.
    Writes out new files to filtered directory
    """
    cand_brcd_indx = candidates_list[1][0]  # candidate barcode indices
    data_obj.m = map(data_obj.m.__getitem__, cand_brcd_indx)
    data_obj.bcs = map(data_obj.bcs.__getitem__, cand_brcd_indx)
    data_obj.save_mex(outfile, compress=True, legacy=legacy)
    print(f"Wrote filtered files to {outfile}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Converts mtx output to h5ad.")

    parser.add_argument("-m", "--mtx_dir", dest="mtx_dir", help="Path to mtx directory.")
    parser.add_argument("-o", "--out_dir", dest="out_dir", help="Path to output directory.")
    ### Optional arguments
    parser.add_argument(
        "-v",
        "--version",
        dest="version",
        default=3,
        help="Version for features file, defines expected file name (1=kallisto,<v3 cellranger,3=star,>=v3 cellranger).",
    )
    parser.add_argument(
        "-d", "--umi_frac_median", dest="umi_frac_median", default=0.01, help="Minimum UMI Fraction of median."
    )
    parser.add_argument("-n", "--umis_nonambient", dest="umis_nonambient", default=500, help="Minimum UMIs nonambient.")
    parser.add_argument("-p", "--max_pvalue", dest="max_pvalue", default=0.01, help="Maximum adjusted p-value.")
    parser.add_argument(
        "-a", "--aligner", dest="aligner", default="cellranger", help="Aligner used: cellranger,alevin,star,kallisto"
    )

    args = vars(parser.parse_args())
    filter_args = [float(args["umi_frac_median"]), float(args["umis_nonambient"]), float(args["max_pvalue"])]
    matrix = None
    legacy = False
    ### Read in all files based on "version", expectation of features/genes names file
    if int(args["version"]) != 3:
        ### kallisto runs here
        ### <v3 cellranger runs here
        matrix = CountMatrix.from_legacy_mtx(args["mtx_dir"])
        legacy = True
    else:
        ### star runs here
        ### >=v3 cellranger runs here
        matrix = CountMatrix.from_v3_mtx(args["mtx_dir"])
    outfile = os.path.join(args["mtx_dir"], "edfiltered")
    matrix.save_mex(outfile, compress=True, legacy=legacy)
    ### Identify candidate barcodes
    candidate_list = find_candidates(matrix, filter_params=filter_args)
    ### Perform Filtering of matrix and write out filtered mtx
    if candidate_list != None:
        filter_candidates(matrix, candidate_list, args["out_dir"], legacy)
    else:
        print("No filtering performed, candidate list is None")
