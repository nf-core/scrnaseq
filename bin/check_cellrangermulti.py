#!/usr/bin/env python

import csv
import os
import sys

def parse_samplesheet(samplesheet_path):
    # Define required headers
    required_headers = ["sample", "multiplexed_sample_id", "description"]

    # Define output directories
    cmo_output_dir = "cmo_files"
    frna_output_dir = "frna_files"
    os.makedirs(cmo_output_dir, exist_ok=True)
    os.makedirs(frna_output_dir, exist_ok=True)

    with open(samplesheet_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        headers = reader.fieldnames

        # Check if required headers are present
        for header in required_headers:
            if header not in headers:
                print(f"Error: Required header '{header}' not found in the samplesheet.")
                return False

        # Process each row
        for row in reader:
            sample = row["sample"]
            multiplexed_sample_id = row["multiplexed_sample_id"]
            description = row["description"]

            # Process CMOs
            if "cmo_ids" in headers and row["cmo_ids"]:
                cmo_filename = os.path.join(cmo_output_dir, f"{sample}_cmo.csv")
                with open(cmo_filename, 'a', newline='') as cmo_file:
                    cmo_writer = csv.writer(cmo_file)
                    if not os.path.exists(cmo_filename) or os.stat(cmo_filename).st_size == 0:
                        cmo_writer.writerow(["sample_id", "cmo_ids", "description"])
                    cmo_writer.writerow([multiplexed_sample_id, row["cmo_ids"], description])

            # Process FRNAs
            if "probe_barcode_ids" in headers and row["probe_barcode_ids"]:
                frna_filename = os.path.join(frna_output_dir, f"{sample}_frna.csv")
                with open(frna_filename, 'a', newline='') as frna_file:
                    frna_writer = csv.writer(frna_file)
                    if not os.path.exists(frna_filename) or os.stat(frna_filename).st_size == 0:
                        frna_writer.writerow(["sample_id", "probe_barcode_ids", "description"])
                    frna_writer.writerow([multiplexed_sample_id, row["probe_barcode_ids"], description])

    print("Parsing completed successfully.")
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_samplesheet>")
        sys.exit(1)

    samplesheet_path = sys.argv[1]
    if not os.path.isfile(samplesheet_path):
        print("Error: Samplesheet file not found.")
        sys.exit(1)

    if not parse_samplesheet(samplesheet_path):
        sys.exit(1)
