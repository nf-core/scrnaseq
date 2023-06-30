#!/usr/bin/env python
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate the lib.csv for cellranger-arc.")

    parser.add_argument("-f", "--fasta", dest="fasta", help="Name of the fasta file.")
    parser.add_argument("-g", "--gtf", dest="gtf", help="Name of the gtf file.")
    parser.add_argument("-m", "--motifs", dest="motifs", help="Name of the motifs file.")
    parser.add_argument("-a", "--add", dest="add", help="Additional filter line.")

    args = vars(parser.parse_args())

    print(args)

    config = open("config", "w")
    config.write("{\n")
    config.write("\torganism: \"{}\"\n".format(args["fasta"].split(".")[0]))
    config.write("\tgenome: [\"cellrangerarc_reference\"]\n")
    config.write("\tinput_fasta: [\"{}\"]\n".format(args["fasta"]))
    config.write("\tinput_gtf: [\"{}\"]\n".format(args["gtf"]))
    config.write("\tinput_motifs: \"{}\"\n".format(args["motifs"]))
    if(args["add"] != "none"):
        config.write(args["add"] + "\n")
    config.write("}")
    config.close()

    print("Wrote config file")
