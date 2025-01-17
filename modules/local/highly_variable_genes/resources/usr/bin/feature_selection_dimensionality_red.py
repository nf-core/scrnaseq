#!/usr/bin/env python3
# ====================================================================================================================
#                                          PRELIMINARIES
# ====================================================================================================================

# MODULE IMPORT

import warnings
warnings.filterwarnings("ignore")


import argparse                     # command line arguments parser
import os                           # filesystem utilities
import scanpy as sc                 # single-cell data processing
import anndata as ad                # store annotated matrix as anndata object
import matplotlib.pyplot as plt     # library for visualization
import seaborn as sns               # library for statistical data visualization
import pandas as pd                 # library for data analysis and manipulation 
import pathlib                      # library for handle filesystem paths


# PARAMETERS

# set script version number
VERSION = "0.0.1"


# ====================================================================================================================
#                                          MAIN FUNCTION
# ====================================================================================================================

def main():
    """
    This function perfors dimensionality reduction of dataset.
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

    parser = argparse.ArgumentParser(prog='DimRed', usage='%(prog)s [options]',description = "Feature selection and dimensionality reduction",
        epilog = "This function reduce the dimensionality of the dataset and only include the most informative genes.")
    parser.add_argument('-ad','--input-h5ad-file',metavar= 'H5AD_INPUT_FILES', type=pathlib.Path, dest='input_h5ad_files',
                        required=True, help="paths of existing count matrix files in h5 format (including file names)")
    parser.add_argument('-o', '--out', metavar='H5AD_OUTPUT_FILE', type=pathlib.Path, default="combined_matrix_DR.h5ad",
                        help="name of the output h5ad file after dimensionality reduction")
    parser.add_argument('-csv', '--csv_out', metavar='CSV_TABLE',type=pathlib.Path, default="UMAP_coordinates.csv",
                        help="csv tabel with UMAP coordinates for each cell")
    parser.add_argument('-r','--results', type=pathlib.Path, default=pathlib.Path('./'), 
                        help="directory to save the results files (default is the current directory)")
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    args = parser.parse_args()

# --------------------------------------------------------------------------------------------------------------------
#                                 DEFINE SAMPLES AND MTX PATHS
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== INPUT H5AD FILES =====")
    input_h5ad_file = args.input_h5ad_files
    output = args.out
    output_csv=args.csv_out

    
    # print info on the available matrices
    print("Reading combined count matrix from the following file:")
    print("-File {}:".format(str(input_h5ad_file)))

# --------------------------------------------------------------------------------------------------------------------
#                                 READ H5AD FILES 
# --------------------------------------------------------------------------------------------------------------------

    
     # Read folders with the MTX combined count matrice and store datasets in a dictionary
    
    print("\n===== READING COMBINED MATRIX =====")
    # read the count matrix for the combined samples and print some initial info
    print("\nProcessing count matrix in folder ... ", end ='')

    adata= sc.read_h5ad(input_h5ad_file)
                
    print("Done!")
    print("Count matrix for combined samples has {} cells and {} genes".format(adata.shape[0],adata.shape[1]))

# --------------------------------------------------------------------------------------------------------------------
#                                 FEATURE SELECTION & DIMENSIONALITY REDUCTION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== FEATURE SELECTION =====")
    # select highly=variable genes for each sample
    print("\nSelecting highly-variable genes are selected within each batch separately and merged")

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'run_id',subset=False)

    print("\n===== DIMENSIONALITY REDUCTION =====")
    print("\nPerforming dimensionality reduction by running principal component analysis (PCA)")
    sc.tl.pca(adata,use_highly_variable = True, n_comps=50)

# --------------------------------------------------------------------------------------------------------------------
#                                 DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION
# --------------------------------------------------------------------------------------------------------------------

    print("\n===== NEAREST NEIGHBOR GRAPH CONSTRUCTION =====")
    print("\nConstruction of the nearest neighbor graph")
    sc.pp.neighbors(adata,n_neighbors=15)

    print("\n===== DIMENSIONALITY REDUCTION FOR DATA VISUALIZATION=====")
    print("\nPerforming dimensionality reduction by running uniform manifold approximation and projection (UMAP)")
    sc.tl.umap(adata,min_dist=0.5)
    
# --------------------------------------------------------------------------------------------------------------------
#                           VISUALIZE UMAP PLOT
# --------------------------------------------------------------------------------------------------------------------

    # Visualize UMAP plot 

    print("\nVisualized UMAP plot")
    sc.pl.umap(adata, color ='run_id',legend_loc='on data',show=False)
    plt.savefig(os.path.join(args.results,'UMAP_plot.png'))
    plt.close()

# --------------------------------------------------------------------------------------------------------------------
#                           SAVE OUTPUT FILE
# --------------------------------------------------------------------------------------------------------------------


    print("\n===== SAVING OUTPUT FILE =====")

    print("Saving h5ad data to file {}".format(output))
    adata.write(output)
    print("Done!")

    df = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names).rename(columns={0: "X_UMAP", 1: "Y_UMAP"})
    df.index.name = 'cell_barcodes'
    print("Saving csv table with UMAP coordinates for each cell {}".format(output_csv))
    df.to_csv(output_csv)
    print("Done!")


#####################################################################################################


if __name__ == '__main__':
    main()