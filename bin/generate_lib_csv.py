#!/usr/bin/env python
import argparse
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate the lib.csv for cellranger-arc.")

    parser.add_argument("-t", "--sample_types", dest="sample_types", help="Comma seperated list of sample types.")
    parser.add_argument("-n", "--sample_names", dest="sample_names", help="Comma seperated list of sample names.")
    parser.add_argument("-f", "--fastq_folder", dest="fastq_folder", help="Folder of FASTQ files.")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")

    args = vars(parser.parse_args())

    print(args)

    sample_types = args["sample_types"].split(",")
    sample_names = args["sample_names"].split(",")
    unique_samples_names = set(sample_names)

    lib_csv = open(args["out"], "w")
    lib_csv.write("fastqs,sample,library_type")

    for i in range(0, len(sample_types)):
        if sample_names[i] in unique_samples_names:
            unique_samples_names.remove(
                sample_names[i]
            )  # this has to be done to account for different Lane files (e.g., L002)
            if sample_types[i] == "gex":
                lib_csv.write("\n{},{},{}".format(args["fastq_folder"], sample_names[i], "Gene Expression"))
            else:
                lib_csv.write("\n{},{},{}".format(args["fastq_folder"], sample_names[i], "Chromatin Accessibility"))

    lib_csv.close()

    print("Wrote lib.csv file to {}".format(args["out"]))
