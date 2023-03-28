import os
import re
import csv
import argparse

parser = argparse.ArgumentParser(description='Convert fastq file names to symlink format')
parser.add_argument('input_file', help='input csv file containing fastq file paths')
parser.add_argument('output_file', help='output csv file containing updated symlink file paths')
args = parser.parse_args()

def create_symlink(fastq_1, fastq_2, sample, lane_num):
    symlink_1 = f"{sample}_S1_L{lane_num}_R1_001.fastq.gz"
    symlink_2 = f"{sample}_S1_L{lane_num}_R2_001.fastq.gz"

    os.symlink(fastq_1, symlink_1)
    os.symlink(fastq_2, symlink_2)

    return symlink_1, symlink_2

def update_samplesheet(samplesheet_file, output_file):
    with open(samplesheet_file, 'r') as f, open(output_file, 'w', newline='') as o:
        reader = csv.DictReader(f)
        writer = csv.DictWriter(o, fieldnames=reader.fieldnames)
        writer.writeheader()

        samples = {}
        for row in reader:
            sample = row['sample']
            fastq_1 = row['fastq_1']
            fastq_2 = row['fastq_2']

            lane_num = 1
            if sample in samples:
                lane_num = samples[sample] + 1
            samples[sample] = lane_num

            symlink_1, symlink_2 = create_symlink(fastq_1, fastq_2, sample, lane_num)
            row['fastq_1'] = os.path.abspath(symlink_1)
            row['fastq_2'] = os.path.abspath(symlink_2)
            writer.writerow(row)

if __name__ == '__main__':
    update_samplesheet(args.input_file, args.output_file)