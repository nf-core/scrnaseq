#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import pathlib                      # library for handle filesystem paths
import scanpy as sc                 # single-cell data processing
from mudata import MuData


warnings.filterwarnings("ignore")

# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function creates a MuData object.
    """
# --------------------------------------------------------------------------------------------------------------------
#                                          LIBRARY CONFIG
# --------------------------------------------------------------------------------------------------------------------

    sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

# --------------------------------------------------------------------------------------------------------------------
#                                          INPUT FROM COMMAND LINE
# --------------------------------------------------------------------------------------------------------------------

# Define command line arguments with argparse

    parser = argparse.ArgumentParser(prog='Create MuData object', usage='%(prog)s [options]',description = "MuData object convertion",
                        epilog = "This function creates a MuData object for storing GEX,VDJ and CITE-seq data.")
    parser.add_argument('-ad','--input-gex-file',metavar= 'GEX_INPUT_FILES', type=pathlib.Path, dest='input_gex_files',
                        help="paths of existing count matrix files in h5ad format (including file names)")
    parser.add_argument('-ai', '--input-vdj-file', metavar='VDJ_INPUT_FILES',type=pathlib.Path, dest='input_vdj_files',
                        default=pathlib.Path(''),help="paths of existing vdj matrix files in h5ad format (including file names)")
    parser.add_argument('-o', '--out', metavar='MUDATA_OUTPUT_FILE', type=pathlib.Path, default="matrix.mudata.h5mu",
                        help="name of the muData object")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT GEX and VDJ FILES =====")
    input_gex_file = args.input_gex_files
    input_vdj_file = args.input_vdj_files
    output = args.out

    # print info on the available matrices
    print("Reading combined gex count matrix from the following file:")
    print(f"-File {input_gex_file}")

    print("Reading filtered annotation table from the following file:")
    print(f"-File {input_vdj_file}")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ GEX AND AB FILES
# --------------------------------------------------------------------------------------------------------------------
    if input_gex_file:
        # Read folders with the MTX combined count matrice and store datasets in a dictionary
        print("\n===== READING COMBINED MATRIX =====")
        # read the gex count matrix for the combined samples and print some initial info
        print("\nProcessing count matrix in folder ... ", end ='')
        adata= sc.read_h5ad(input_gex_file)
        print("Done!")
        print(f"Gex count matrix for combined samples has {adata.shape[0]} cells and {adata.shape[1]} genes")
    else:
        print("No valid input file provided. Skipping reading of the count matrix.")

# --------------------------------------------------------------------------------------------------------------------
#                                 READ VDJ FILES
# --------------------------------------------------------------------------------------------------------------------
    if input_vdj_file and input_vdj_file != pathlib.Path(''):
        # Read folders with the filtered contigue annotation and store datasets in a dictionary
        print("\n===== READING CONTIGUE ANNOTATION MATRIX =====")
        # read the filtered contigue annotation file for the combined samples and print some initial info
        print("\nProcessing filtered contigue table in folder ... ", end ='')
        adata_vdj= sc.read_h5ad(input_vdj_file)
        print("Done!")
    else:
        print("No valid input file provided. Skipping reading of the vdj annotation.")

# --------------------------------------------------------------------------------------------------------------------
#                                 CREATE MUDATA OBJECT
# --------------------------------------------------------------------------------------------------------------------
    #Creates dictionary to store all modalities
    modalities = {}
    try:
        # Add 'gex' modality if defined
        if adata[:, adata.var["feature_types"] == "Gene Expression"].shape[1] > 0:
            modalities["gex"] = adata[:, adata.var["feature_types"] == "Gene Expression"]
        # Add 'pro' modality if defined
        if adata[:, adata.var["feature_types"] == "Antibody Capture"].shape[1] > 0:
            modalities["pro"] = adata[:, adata.var["feature_types"] == "Antibody Capture"]
    except NameError:
        pass

    try:
        # Add 'airr' modality if defined
        if adata_vdj is not None:
            modalities["airr"] = adata_vdj
    except NameError:
        pass

    # Creates MuData object
    mdata = MuData(modalities)

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== SAVING OUTPUT FILE =====")

    print(f"Saving MuData object to file {output}")
    mdata.write(output)
    print("Done!")

#####################################################################################################

if __name__ == '__main__':
    main()
