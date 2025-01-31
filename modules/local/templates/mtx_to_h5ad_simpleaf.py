#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
import argparse
import anndata
from anndata import AnnData
import platform
import json

def say(quiet, words):
    if not quiet:
        print(words)

def load_fry(
    frydir,
    output_format="snRNA",
    aux_columns=["X", "Y"],
    gene_id_to_name=None,
    quiet=False,
):
    """
    load alevin-fry quantification result into an AnnData object

    Required Parameters
    ----------
    frydir : `str`
        The path to a output directory returned by alevin-fry quant command. \\
        The directory containing the alevin-fry quantification (i.e. the the quant.json file & alevin subdirectory).

    Optional Parameters
    ----------
    output_format : `str` or `dict`
        A string represents one of the pre-defined output formats, which are "scRNA", "S+A", "snRNA", "all", "U+S+A" and "velocity". \\
        If a customized format of the returned `AnnData` is needed, one can pass a Dictionary.\\
        See Notes section for details.

    aux_columns : `list[str]` (default: `["X", "Y"]`)
        A list of strings contains the column names of the auxiliary information in the barcodes file starting from the second column. \\
        The first column is assumed to be the barcodes and is named as "barcodes". \\
        Extra auxiliary columns in the barcodes file without a specified name will be ignored.

    gene_id_to_name : `str` or `None` (default: `None`)
        The path to a file that contains the mapping from gene names to gene ids. \\
        It is only needed if \\
            1. you are not using the simpleaf pipeline (`simpleaf index` + `simpleaf quant`), \\
            2. you have such a file, and,
            3. you want to add this information to the coldata of your anndata.
        If you do, please ensure it is a tab-separated, two-column file without a header, and the first column is the gene ids and the second column is the gene names.

    nonzero : `bool` (default: `False`)
        True if cells with non-zero expression value across all genes should be filtered in each layer.
        False if unexpressed genes should be kept.

    quiet : `bool` (default: `False`)
        True if function should be quiet.
        False if messages (including error messages) should be printed out.

    Notes
    ----------
    The `output_format` argument takes either a dictionary that defines the customized format or
    a string that represents one of the pre-defined format of the returned `AnnData` object.

    Each of the pre-defined formats contains a `X` field and some optional extra `AnnData.layers`
    obtained from the submatrices representing unspliced (U), spliced (S) and ambiguous (A) counts
    returned by alevin-fry.

    The following formats are defined:

    * "scRNA" and "S+A": \\
        This format is recommended for single cell RNA-sequencing experiments.
        It returns a `X` field that contains the S+A count of each gene in each cell,
        and a `unspliced` field that contains the U count of each gene.

    * "snRNA", "all" and "U+S+A": \\
        These three formats are the same. They return a `X` field that contains the U+S+A
        count of each gene in each cell without any extra layers.
        It is recommended for single-nucleus RNA-sequencing experiments.
        CellRanger 7 returns this format for both single-cell and single-nucleus experiments.

    * "raw": \\
        This format uses the S count matrix as the `X` field and put the U, S, and A counts into three
        separate layers, which are "unspliced", "spliced" and "ambiguous".

    * "velocity": \\
        This format is the same as "scRNA", except it contains two extra layers: the "spliced" layer,
        which contains the S+A counts, and the "unspliced" layer, which contains the U counts.

    A custom output format can be defined using a Dictionary specifying the desired format of the output `Anndata` object.
    If the input is not a USA mode quantification directory, this parameter is ignored
    and the count matrix is returned in the `X` field of the returned `AnnData` object.  If the input
    quantification directory contains a USA mode quantification, then there are 3 sub-matrices that can
    be referenced in the dictionary; 'U', 'S', 'A' containing, respectively, unspliced, spliced and
    ambiguous counts.  The dictionary should have entries of the form `key` (str) : `value` (list[str]).
    The following constraints apply : there should be one key-value pair with the key `X`, the resulting
    counts will be returned in the `X` field of the AnnData object. There can be an arbitrary number
    of other key-value pairs, but each will be returned as a layer of the resulting AnnData object.
    Within the key-value pairs, the key refers to the layer name that will be given to the combined
    count matrix upon output, and the value should be a subset of `['U', 'S', 'A']` that defines
    which sub-matrices should be summed.  For example:
    `{'X' : ['S', 'A'], 'unspliced' : ['U']}`
    will result in a return AnnData object where the X field has a matrix in which each entry
    corresponds to the summed spliced and ambiguous counts for each gene in each cell, and there
    is an additional "unspliced" layer, whose counts are taken directly from the unspliced sub-matrix.

    Returns:
    ----------
        An AnnData object with X and layers corresponding to the requested `output_format`.

    """

    # since alevin-fry 0.4.1 the generic "meta_info.json"
    # has been replaced by a more informative name for each
    # sub-command. For quantification, it is "quant.json".
    # we check for both files here, in order.
    meta_info_files = ["quant.json", "meta_info.json"]

    fpath = os.path.sep.join([frydir, meta_info_files[0]])
    # first, check for the new file, if we don't find it, check
    # for the old one.
    if not os.path.exists(fpath):
        say(
            quiet,
            f"Did not find a {meta_info_files[0]} file, checking for older {meta_info_files[1]}.",
        )

        fpath = os.path.sep.join([frydir, meta_info_files[1]])
        # if we don't find the old one either, then return None
        if not os.path.exists(fpath):
            raise IOError(
                "The profvided `frydir` doesn't contain required meta info file; cannot proceed."
            )

    # if we got here then we had a valid json file, so
    # use it to get the number of genes, and if we are
    # in USA mode or not.
    meta_info = json.load(open(fpath))
    ng = meta_info["num_genes"]
    usa_mode = meta_info["usa_mode"]

    say(quiet, f"USA mode: {usa_mode}")

    # if we are in USA mode
    if usa_mode:
        # preparation
        # each gene has 3 splicing statuses, so the actual number of distinct
        # genes is ng/3.
        ng = int(ng / 3)
        output_assays = process_output_format(output_format, quiet)
    else:
        say(
            quiet,
            "Processing input in standard mode, the count matrix will be stored in field 'X'.",
        )
        if output_format != "scRNA":
            say(quiet, "Output_format will be ignored.")

    # read the gene ids
    afg_df = pd.read_table(
        os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"]),
        names=["gene_id"],
        nrows=ng
    )

    # if we have a gene name to id mapping, use it
    # Otherwise, we see if there is a default mapping file
    # We first hope the user gives one
    gene_id_to_name_path = gene_id_to_name
    # make sure the file exists
    if gene_id_to_name_path is not None:
        if not os.path.exists(gene_id_to_name_path):
            say(
                quiet,
                f"The provided gene id to name mapping file {gene_id_to_name_path} does not exist; ignored.",
            )
            gene_id_to_name_path = None

    # we then check if we can find a default mapping file
    if gene_id_to_name_path is None:
        default_gene_id_to_name_path = os.path.sep.join([frydir, "gene_id_to_name.tsv"])
        if os.path.exists(default_gene_id_to_name_path):
            gene_id_to_name_path = default_gene_id_to_name_path
            say(
                quiet,
                f"Using simpleaf gene id to name mapping file: {gene_id_to_name_path}",
            )

    # read the file if we find it
    if gene_id_to_name_path is not None:
        gene_id_to_name_df = pd.read_table(
            gene_id_to_name_path, names=["gene_id", "gene_symbol"]
        )
        afg_df = pd.merge(afg_df, gene_id_to_name_df, on="gene_id", how="left")

    afg_df = afg_df.set_index("gene_id", drop=False)

    # read the barcodes file
    abc_df = pd.read_table(
        os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"]), header=None
    )

    # if column names and num columns don't match, match them
    # we have only barcode
    if abc_df.shape[1] == 1:
        columns = ["barcodes"]
    else:
        say(quiet, f"Found {abc_df.shape[1] - 1} auxiliary columns in barcodes file.")
        columns = ["barcodes"] + aux_columns

    if abc_df.shape[1] != len(columns):
        ncol = min(abc_df.shape[1], len(columns))

        say(
            quiet,
            f"Number of auxiliary columns in barcodes file does not match the provided column names. Using the first {ncol} columns in the barcode file with names: {columns[:ncol]}.",
        )

        abc_df.columns = columns[:ncol]
        columns = columns[:ncol]

    abc_df = abc_df.set_axis(columns, axis=1)
    abc_df = abc_df.set_index(columns[0], drop=False)

    say(quiet, "Reading the quantification matrix.")
    af_raw = sc.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
    x = af_raw.X
    # if we're not in USA mode, just combine this info into
    # an AnnData object
    if not usa_mode:
        say(quiet, "Constructing the AnnData object.")
        af = sc.AnnData(x.T, var=abc_df, obs=afg_df)
        af = af.T

    else:  # USA mode
        say(quiet, "Converting USA mode quantification into specified layers.")
        # otherwise, combine the sub-matrices into the output object as
        # specified by `output_assays`
        rd = {"S": range(0, ng), "U": range(ng, 2 * ng), "A": range(2 * ng, 3 * ng)}
        xcounts = output_assays["X"]
        o = x[:, rd[xcounts[0]]]
        for wc in xcounts[1:]:
            o += x[:, rd[wc]]
        af = sc.AnnData(o.T, var=abc_df, obs=afg_df)
        af = af.T

        # now, if there are other layers requested, populate those
        for other_layer in output_assays.keys() - "X":
            xcounts = output_assays[other_layer]
            o = x[:, rd[xcounts[0]]]
            for wc in xcounts[1:]:
                o += x[:, rd[wc]]
            af.layers[other_layer] = o

    say(quiet, "Done.")
    return af


def process_output_format(output_format, quiet):
    # make sure output_format isn't empty
    if not output_format:
        raise ValueError("output_format cannot be empty")

    if isinstance(output_format, (str, dict)):
        if isinstance(output_format, str):
            predefined_format = {
                "scrna": {"X": ["S", "A"], "unspliced": ["U"]},
                "S+A": {"X": ["S", "A"]},
                "snrna": {"X": ["U", "S", "A"]},
                "all": {"X": ["U", "S", "A"]},
                "U+S+A": {"X": ["U", "S", "A"]},
                "velocity": {
                    "X": ["S", "A"],
                    "spliced": ["S", "A"],
                    "unspliced": ["U"],
                },
                "raw": {
                    "X": ["S"],
                    "spliced": ["S"],
                    "unspliced": ["U"],
                    "ambiguous": ["A"],
                },
            }

            output_format = output_format.lower()
            if output_format not in predefined_format.keys():
                # invalid output_format string
                say(quiet, "A undefined Provided output_format string provided.")
                say(quiet, "See function help message for details.")
                raise ValueError("Invalid output_format.")
            say(quiet, f"Using pre-defined output format: {output_format}")
            say(
                quiet,
                f"Will populate output field X with sum of counts from {predefined_format[output_format]['X']}.",
            )
            for k, v in predefined_format[output_format].items():
                if k != "X":
                    say(quiet, f"Will combine {v} into output layer {k}.")

            return predefined_format[output_format]
        else:
            say(quiet, "Processing user-defined output format.")
            # make sure the X is there
            if "X" not in output_format.keys():
                raise ValueError(
                    'In USA mode some sub-matrices must be assigned to the "X" (default) output.'
                )
            say(
                quiet,
                f"Will populate output field X with sum of counts frorm {output_format['X']}.",
            )

            for k, v in output_format.items():
                if not v:
                    # empty list
                    raise ValueError(
                        f"The element list of key '{k}' in output_format is empty. Please remove it."
                    )

                # v contains Non-USA element
                if len(set(v) - set(["U", "S", "A"])) != 0:
                    # invalid value
                    raise ValueError(
                        f"Found non-USA element in output_format element list '{v}' for key '{k}'; cannot proceed."
                    )
                if k != "X":
                    say(quiet, f"Will combine {v} into output layer {k}.")

            return output_format
    else:
        raise ValueError(
            "Provided invalid output_format. See function help message for details"
        )

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.
    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.
    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

def dump_versions():
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "pandas": pd.__version__,
            "anndata": anndata.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))

def input_to_adata(
    input_data: str,
    output: str,
    sample: str,
):
    print(f"Reading in {input_data}")

    # open main data
    simpleaf_h5ad_path = f"{input_data}/alevin/quant.h5ad"

    if os.path.exists(simpleaf_h5ad_path):
        adata = sc.read_h5ad(simpleaf_h5ad_path)
    else:
        output_format={
            "X": ["S", "U", "A"],
            "spliced": ["S"],
            "unspliced": ["U"],
            "ambiguous": ["A"],
        }
        adata = load_fry(input_data, output_format=output_format, quiet=False)

    adata.obs["sample"] = sample

    # standard format
    # index are gene IDs and symbols are a column
    adata.var['gene_versions'] = adata.var.index
    adata.var.index = adata.var['gene_versions'].str.split('.').str[0].values
    adata.var_names_make_unique()

    # write results
    adata.write_h5ad(f"{output}")
    print(f"Wrote h5ad file to {output}")

#
# Run main script
#

# create the directory with the sample name
os.makedirs("${meta.id}", exist_ok=True)

# input_type comes from NF module
input_to_adata(
    input_data="${inputs}",
    output="${meta.id}_${meta.input_type}_matrix.h5ad",
    sample="${meta.id}"
)

# dump versions
dump_versions()
