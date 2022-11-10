#!/usr/bin/env python
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate the lib.csv for cellranger-arc.")

    parser.add_argument("-t", "--sample_types", dest="in", help="Comma seperated list of sample types.")
    parser.add_argument("-n", "--samples_names", dest="in", help="Comma seperated list of sample names.")
    parser.add_argument("-o", "--out", dest="out", help="Output path.")

    args = vars(parser.parse_args())

    print(args)

    sample_types = args["sample_types"].split(",")
    sample_names = args["samples_names"].split(",")
    sample = args["sample"]

    lib_csv = open(args["out"], "w")
    lib_csv.write("fastqs,sample,library_type")


    for i in range(0,len(sample_types)):
        if(sample_types[i] == "gex"):
            lib_csv.write(".,{},{}".format(sample_names[i],"Gene Expression"))
        else:
            lib_csv.write(".,{},{}".format(sample_names[i],"Chromatin Accessibility"))

    lib_csv.close()

    print("Wrote lib.csv file to {}".format(args["out"]))
