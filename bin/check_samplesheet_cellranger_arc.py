#!/usr/bin/env python

"""Provide a command line tool to validate and transform tabular samplesheets for multiome data."""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical("The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure`_::

            sample,folder_ATAC,folder_GEX,expected_cells,seq_center
            SAMPLE_1,SAMPLE_1_ATAC,SAMPLE_2_GEX,n,name
            SAMPLE_2,SAMPLE_1_ATAC,SAMPLE_2_GEX,n,name
            SAMPLE_3,SAMPLE_1_ATAC,SAMPLE_2_GEX,n,name

    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 2
        MIN_HEADER = ["sample", "folder_ATAC", "folder_GEX"]
        OPT_HEADER = ["expected_cells", "seq_center"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        unknown_header = 0
        min_header_count = 0
        colmap = {"sample": 0, "folder_ATAC": 1, "folder_GEX": 2}
        i = 0
        for h in header:
            if h not in MIN_HEADER and h not in OPT_HEADER:
                unknown_header = 1
            if h in MIN_HEADER:
                min_header_count = min_header_count + 1
            colmap[h] = i
            i = i + 1
        if unknown_header or min_header_count < len(MIN_HEADER):
            given = ",".join(header)
            wanted = ",".join(HEADER)
            print(f"ERROR: Please check samplesheet header -> {given} != {wanted}")
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(header):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(header)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, folder_ATAC, folder_GEX = lspl[: len(MIN_HEADER)]
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check expected cells is an integer if present
            expected_cells = ""
            if "expected_cells" in header:
                expected_cells = lspl[colmap["expected_cells"]]
                if not is_integer(expected_cells):
                    print_error("Expected cells must be an integer", "Line", line)

            ## If present, replace spaces with _ in sequencing center name
            seq_center = ""
            if "seq_center" in header:
                seq_center = lspl[colmap["seq_center"]]
                seq_center = seq_center.replace(" ", "_")

            ## Create sample mapping dictionary = { sample: [ folder_ATAC, folder_GEX ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    # print_error("Samplesheet contains duplicate rows!", "Line", line)
                    sample_mapping_dict[sample].append(sample_info)
                else:
                    sample_mapping_dict[sample].append(sample_info)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
