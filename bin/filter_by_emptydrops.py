#!/usr/bin/env python
import argparse
import sys
import os
from typing import NamedTuple
from xmlrpc.client import Boolean
import pandas as pd
import scanpy as sc
import scipy.sparse as sp_sparse
import scipy.io as sp_io
from emptydrops import find_nonambient_barcodes
from emptydrops.matrix import CountMatrix
import emptydrops.h5_constants as h5_constants
from emptydrops.feature_ref import FeatureReference, FeatureDef
import biomart  


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
        """
        # from vpolo.alevin import parser
        # alevin https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py
        
        barcodes_tsv = os.path.join(genome_dir, "quants_mat_rows.txt")
        features_tsv = os.path.join(genome_dir, "quants_mat_cols.txt")
        matrix_mtx = os.path.join(genome_dir, "quants_mat.mtx")
        for filepath in [barcodes_tsv, features_tsv, matrix_mtx]:
            if not os.path.exists(filepath):
                raise IOError("Required file not found: %s" % filepath)
        barcodes = pd.read_csv(barcodes_tsv, delimiter='\t', header=None, usecols=[0]).values.squeeze()
        genes = pd.read_csv(features_tsv, delimiter='\t', header=None, usecols=[0], names=['gene_id'])
        # list_ids = [row['gene_id'] for idx, row in genes.iterrows()]
        # gene_map = get_ensembl_mappings(list_ids)
        feature_defs = [FeatureDef(idx, row['gene_id'], row['gene_id'], "Gene Expression", []) for idx, row in genes.iterrows()]

        # get_ensembl_mappings(feature_defs)
        feature_ref = FeatureReference(feature_defs, [])
        matrix = sp_io.mmread(matrix_mtx)
        mat = CountMatrix(feature_ref, barcodes, matrix)
        mat.tocsc()
        return mat


    @staticmethod
    def from_kallisto(genome_dir):
        """
        Override original implementation for different aligner inputs
        """
        # how to get the full file name?
        barcodes_tsv = os.path.join(genome_dir, "cells_x_genes.barcodes.txt")
        features_tsv = os.path.join(genome_dir, "cells_x_genes.genes.txt")
        matrix_mtx = os.path.join(genome_dir, "cells_x_genes.mtx")         # designate name of matrix file
        for filepath in [barcodes_tsv, features_tsv, matrix_mtx]:
            if not os.path.exists(filepath):
                raise IOError("Required file not found: %s" % filepath)
        barcodes = pd.read_csv(barcodes_tsv, delimiter='\t', header=None, usecols=[0]).values.squeeze()
        genes = pd.read_csv(features_tsv, delimiter='\t', header=None, usecols=[0], names=['gene_id'])
        feature_defs = [FeatureDef(idx, row['gene_id'], row['gene_id'], "Gene Expression", []) for idx, row in genes.iterrows()]
        feature_ref = FeatureReference(feature_defs, [])
        matrix = sp_io.mmread(matrix_mtx).tocsc()
        mat = CountMatrix(feature_ref, barcodes, matrix)
        return mat

    @staticmethod
    def from_anndata(adata):
        """
        Override original implementation for different aligner inputs
        """
        barcodes = adata.obs_names.values
        genes = adata.var_names.values
        feature_defs = [FeatureDef(idx, gene_id, None, "Gene Expression", []) for (idx, gene_id) in enumerate(genes)]
        feature_ref = FeatureReference(feature_defs, [])
        matrix = adata.X.T.astype(int)
        if type(matrix) is not sp_sparse.csc_matrix:
            matrix = matrix.tocsc()
        mat = CountMatrix(feature_ref, barcodes, matrix)
        return mat

    @staticmethod
    def from_general(genome_dir,aligner):
        """
        Accounts for cellranger and star
        """
        # {aligner:(barcodes,features,matrix)}
        align_map = {
            'cellranger':('barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz',),
            'alevin':('quants_mat_rows.txt','quants_mat_cols.txt','quants_mat.mtx',),
            'star':('barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz',),
            'kallisto':('barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz',),
        }
        barcodes_tsv = os.path.join(genome_dir, "barcodes.tsv.gz")
        features_tsv = os.path.join(genome_dir, "features.tsv.gz")
        matrix_mtx = os.path.join(genome_dir, "matrix.mtx.gz")
        for filepath in [barcodes_tsv, features_tsv, matrix_mtx]:
            if not os.path.exists(filepath):
                raise IOError("Required file not found: %s" % filepath)
        barcodes = pd.read_csv(barcodes_tsv, delimiter='\t', header=None, usecols=[0]).values.squeeze()
        features = pd.read_csv(features_tsv, delimiter='\t', header=None)

        feature_defs = []
        for (idx, (_, r)) in enumerate(features.iterrows()):
            fd = FeatureDef(idx, r[0], r[1], r[2], [])
            feature_defs.append(fd)

        feature_ref = FeatureReference(feature_defs, [])

        matrix = sp_io.mmread(matrix_mtx).tocsc()
        mat = CountMatrix(feature_ref, barcodes, matrix)
        return mat

    @staticmethod
    def from_v3_mtx(genome_dir):
        """
        Accounts for cellranger and star
        """
        barcodes_tsv = os.path.join(genome_dir, "barcodes.tsv.gz")
        features_tsv = os.path.join(genome_dir, "features.tsv.gz")
        matrix_mtx = os.path.join(genome_dir, "matrix.mtx.gz")
        for filepath in [barcodes_tsv, features_tsv, matrix_mtx]:
            if not os.path.exists(filepath):
                raise IOError("Required file not found: %s" % filepath)
        barcodes = pd.read_csv(barcodes_tsv, delimiter='\t', header=None, usecols=[0]).values.squeeze()
        features = pd.read_csv(features_tsv, delimiter='\t', header=None)

        feature_defs = []
        for (idx, (_, r)) in enumerate(features.iterrows()):
            fd = FeatureDef(idx, r[0], r[1], r[2], [])
            feature_defs.append(fd)

        feature_ref = FeatureReference(feature_defs, [])
        matrix = sp_io.mmread(matrix_mtx).tocsc()
        mat = CountMatrix(feature_ref, barcodes, matrix)
        return mat

    @staticmethod
    def load_mtx(mtx_dir):
        """
        Override original implementation for different aligner inputs
        """
        legacy_fn = os.path.join(mtx_dir, "genes.tsv")
        v3_fn = os.path.join(mtx_dir, "features.tsv.gz")

        if os.path.exists(legacy_fn):
            return CountMatrix.from_legacy_mtx(mtx_dir)

        if os.path.exists(v3_fn):
            return CountMatrix.from_v3_mtx(mtx_dir)

        raise IOError("Not a valid path to a feature-barcode mtx directory: '%s'" % str(mtx_dir))


def get_ensembl_mappings(list_ids):                                   
    # Set up connection to server                                               
    server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')         
    # modularize by species?
    mart = server.datasets['mmusculus_gene_ensembl']                            
                                                                                
    # List the types of data we want                                            
    attributes = [
                    'mgi_symbol', 
                    'ensembl_gene_id',
                 ]
                                                            
    # Get the mapping between the attributes
    response = mart.search({'attributes': attributes,
                 'filters': {'ensembl_gene_id': list_ids}})
    data = response.raw.data.decode('ascii')
                                                                                
    ensembl_to_genesymbol = {}                                                                                                    
    for line in data.splitlines():                                              
        line = line.split('\t')                                                 
        # The entries are in the same order as in the `attributes` variable                                                 
        gene_symbol = line[0]                                                   
        ensembl_gene = line[1]                                               
                                                                                
        ensembl_to_genesymbol[ensembl_gene] = gene_symbol                                       
                                                                                
    return ensembl_to_genesymbol

        
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
    try:
        data_obj.save_mex(outfile, compress=True, legacy=legacy)
        print(f"Wrote filtered files to {outfile}")
    except Exception as e:
        print(e)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filters results using emptyDrops.")

    parser.add_argument("-m", "--mtx_dir", dest="mtx_dir", required = True, help="Path to mtx directory.")
    parser.add_argument("-o", "--out_dir", dest="out_dir", required = True, help="Path to output directory.")
    ### Optional arguments
    parser.add_argument(
        "-a",
        "--aligner",
        dest="aligner",
        default= "cellranger",
        help="Aligner used: cellranger,alevin,star,kallisto",
    )    
    parser.add_argument(
        "-d", 
        "--umi_frac_median", 
        dest="umi_frac_median", 
        default=0.01, 
        help="Minimum UMI Fraction of median."
    )
    parser.add_argument(
        "-n", 
        "--umis_nonambient", 
        dest="umis_nonambient", 
        default=500, 
        help="Minimum UMIs nonambient.")
    parser.add_argument(
        "-p", 
        "--max_pvalue", 
        dest="max_pvalue", 
        default=0.01, 
        help="Maximum adjusted p-value.")

    args = vars(parser.parse_args())
    filter_args = [float(args["umi_frac_median"]), float(args["umis_nonambient"]), float(args["max_pvalue"])]
    matrix = None
    legacy = False
    ### Read in all files based on "version", expectation of features/genes names file
    
    if args.get('aligner') in ('kallisto'): ### kallisto runs here
        matrix = GeneralCountMatrix.from_kallisto(args["mtx_dir"])
        legacy = True
    elif args.get('aligner') in ('alevin'): ### alevin runs here
        matrix = GeneralCountMatrix.from_alevin(args["mtx_dir"])
    elif args.get('aligner') in ('cellranger','star'):  ### >=v3 cellranger, star runs here
        matrix = GeneralCountMatrix.from_v3_mtx(args["mtx_dir"])
    
    # #test
    # outfile = os.path.join(args["mtx_dir"], "edfiltered")
    # matrix.save_mex(outfile, compress=True, legacy=legacy)
    ### Identify candidate barcodes
    print(matrix.bcs)
    candidate_list = find_candidates(matrix, filter_params=filter_args)
    ### Perform Filtering of matrix and write out filtered mtx
    if candidate_list != None:
        filter_candidates(matrix, candidate_list, args["out_dir"], legacy)
    else:
        print("No filtering performed, candidate list is None")
