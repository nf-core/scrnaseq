#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT
import warnings
import argparse                     # command line arguments parser
import pathlib                      # library for handle filesystem paths
import glob
import scanpy as sc                 # single-cell data processing
import scirpy as ir                 # single-cell AIRR-data
import anndata as ad                # store annotated matrix as anndata object


warnings.filterwarnings("ignore")

# PARAMETERS
# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function concatenates csv files from vdj modality.
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

    parser = argparse.ArgumentParser(prog='Concatenate_vdj', usage='%(prog)s [options]',description = "VDJ data concatenation",
                        epilog = "This function concatenated vdj filtered contig annotation files into a single csv files.")
    parser.add_argument('-ai', '--input-vdj-dir', metavar='VDJ_INPUT_FILES',nargs='+',type=pathlib.Path, dest='input_vdj_files',
                        help="paths of existing directory containing vdj matrix files in csv format (including file names)")
    parser.add_argument('-id', '--input-run-id', metavar='INPUT_RUN_ID', nargs='+', dest='input_run_id',
                        help="names of the run-id corresponding to the input adata")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="combined.vdj.h5ad",
                        help="name of the h5ad object containing the concatenated vdj table")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== VDJ FILES =====")
    input_vdj_file = args.input_vdj_files
    input_run_id = args.input_run_id
    output = args.out

    # print info on the available matrices
    print("Reading vdj matrix from the following files:")
    for run, mtx in zip(input_run_id, input_vdj_file):
        print(f"Run: {run:15s} - File: {mtx}")


# --------------------------------------------------------------------------------------------------------------------
#                                 READ VDJ FILES
# --------------------------------------------------------------------------------------------------------------------
    
    vdj_files = []
    for folder in glob.glob("*/filtered_contig_annotations.csv"):
        vdj_files.append(folder)

    adata_vdj_list = []
    if vdj_files:

        for run, vdj in zip(input_run_id,vdj_files):
            # Read folders with the filtered contigue annotation and store datasets in a dictionary
            print("\n===== READING CONTIGUE ANNOTATION MATRIX =====")
            print("\nProcessing filtered contigue table in folder ... ", end ='')
            adata_vdj= ir.io.read_10x_vdj(vdj)
            print("Done!")
            adata_vdj_list.append(adata_vdj)
    else:
        print("No valid input file provided. Skipping reading of the vdj annotation.")

# --------------------------------------------------------------------------------------------------------------------
#                           VDJ TABLE CONCATENATION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== CONCATENATING VDJ TABLES =====")

    if len(adata_vdj_list) > 1:
        adata_vdj_concatenated = ad.concat(adata_vdj_list, join= "outer", merge ="same", label="sample",
                                        keys= input_run_id, index_unique="_")

    print(f"Concatenated vdj table for {len(input_run_id)} batched has {adata_vdj_concatenated.shape[0]} cells")
    print("Done!")
# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== SAVING OUTPUT FILE =====")

    print(f"Saving vdj table data in {output}")
    adata_vdj_concatenated.write(output)
    print("Done!")


#####################################################################################################


if __name__ == '__main__':
    main()
