#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import warnings
warnings.filterwarnings("ignore")

import pathlib                      # library for handle filesystem paths
import argparse                     # command line arguments parser
import scanpy as sc                 # single-cell data processing

# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function lognormalized the matrices in h5ad format.
    """
# --------------------------------------------------------------------------------------------------------------------
#                                          LIBRARY CONFIG
# --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

# Define command line arguments with argparse

    parser = argparse.ArgumentParser(prog='LogNorm', usage='%(prog)s [options]',description = "Normalization and logaritmic trasformation of count matrix",
                        epilog = "This function normalize and logarithmize the data")
    parser.add_argument('-ad','--input-h5ad-file',metavar= 'H5AD_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="matrix.norm.h5ad",
                        help="path and name of the output h5ad file after filtering")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()
    parser.print_help()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5ad_files = args.input_h5ad_files
    output = args.out


    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print(f"-File {str(input_h5ad_files)}:")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES
# --------------------------------------------------------------------------------------------------------------------


     # Read folders with the MTX combined count matrice and store datasets in a dictionary

    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end ='')

    adata= sc.read_h5ad(input_h5ad_files)

    print("Done!")
    print(f"Count matrix for combined samples has {adata.shape[0]} cells and {adata.shape[1]} genes")

# --------------------------------------------------------------------------------------------------------------------
#                                 NORMALIZATION
# --------------------------------------------------------------------------------------------------------------------
    #Saving count data before normalization
    print("Saving count data before normalization in slot Count.")
    adata.layers["count"] = adata.X.copy()

    print("\n===== NORMALIZATION =====")
    # Normalizing to median total counts
    print("\nNormalize to median total counts ... ")
    sc.pp.normalize_total(adata)

    print("Done!")

    print("\n===== LOGARITMIC TRASFORMATION =====")
    print("\nLogarithmize the data ... ", end ='')
    # Logarithmize the data
    sc.pp.log1p(adata)

    print("Done!")

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== SAVING OUTPUT FILE =====")

    #Saving count data before normalization
    print("Saving lognormalized data in slot normalized")
    adata.layers["normalized"] = adata.X.copy()
    print("Saving h5ad data to file {}".format(output))
    adata.write(output)
    print("Done!")

#####################################################################################################


if __name__ == '__main__':
    main()

