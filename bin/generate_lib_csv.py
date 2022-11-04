#!/usr/bin/env python
import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate the lib.csv for cellranger-arc.")

    parser.add_argument("-i", "--in", dest="in", help="Path to samplesheet.")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")

    args = vars(parser.parse_args())

    in_file = args["in"]

    print("Wrote lib.csv file to {}".format(args["out"]))
