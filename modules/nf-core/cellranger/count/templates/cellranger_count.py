#!/usr/bin/env python3
"""
Automatically rename staged files for input into cellranger count.

Copyright (c) Gregor Sturm 2023 - MIT License
"""

import re
import shlex
from collections import defaultdict
from pathlib import Path
from subprocess import run
from textwrap import dedent


def chunk_iter(seq, size):
    """iterate over `seq` in chunks of `size`"""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


sample_id = "${meta.id}"

# get fastqs, ordered by path. Files are staged into
#   - "fastq_001/{original_name.fastq.gz}"
#   - "fastq_002/{original_name.fastq.gz}"
#   - ...
fastqs = sorted(Path(".").glob("fastq_*/*"))

if not fastqs:
    raise AssertionError("No FASTQ files found under fastq_*/ for cellranger count.")

# target directory in which the renamed fastqs will be placed
fastq_all = Path("./fastq_all")
fastq_all.mkdir(exist_ok=True)

# Match R1 in the filename, but only if it is followed by a non-digit or non-character
# match "file_R1.fastq.gz", "file.R1_000.fastq.gz", etc. but
# do not match "SRR12345", "file_INFIXR12", etc
filename_pattern = r"([^a-zA-Z0-9])R1([^a-zA-Z0-9])"

# Match SRA-style suffixes such as "_1.fastq.gz", "_2.fq", "_3.fastq"
sra_suffix_pattern = re.compile(r"(_[123])(\.f(ast)?q(\.gz)?)$")

# First, try to group SRA-style *_1/_2/_3 files by prefix
sra_groups = []
seen_in_sra = set()
by_prefix = defaultdict(dict)

for p in fastqs:
    m = sra_suffix_pattern.search(p.name)
    if not m:
        continue
    suffix = m.group(1)  # "_1", "_2" or "_3"
    prefix = p.name[: m.start(1)]
    by_prefix[prefix][suffix] = p

for prefix, files in by_prefix.items():
    # Require at least _2 and _3 to define R1/R2
    if "_2" in files and "_3" in files:
        group = {
            "prefix": prefix,
            "i1": files.get("_1"),
            "r1": files["_2"],
            "r2": files["_3"],
        }
        sra_groups.append(group)
        seen_in_sra.update(files.values())

# Remaining files are treated as bcl2fastq-style R1/R2 pairs
remaining_fastqs = [p for p in fastqs if p not in seen_in_sra]

lane = 1

if remaining_fastqs:
    # For bcl2fastq-style data we still assume [R1, R2, R1, R2, ...]
    assert len(remaining_fastqs) % 2 == 0

    for r1, r2 in chunk_iter(remaining_fastqs, 2):
        # double escapes are required because nextflow processes this python 'template'
        if re.sub(filename_pattern, r"\\1R2\\2", r1.name) != r2.name:
            raise AssertionError(
                dedent(
                    f"""\
                    We expect R1 and R2 of the same sample to have the same filename except for R1/R2.
                    This has been checked by replacing "R1" with "R2" in the first filename and comparing it to the second filename.
                    If you believe this check shouldn't have failed on your filenames, please report an issue on GitHub!

                    Files involved:
                        - {r1}
                        - {r2}
                    """
                )
            )
        r1.rename(fastq_all / f"{sample_id}_S1_L{lane:03d}_R1_001.fastq.gz")
        r2.rename(fastq_all / f"{sample_id}_S1_L{lane:03d}_R2_001.fastq.gz")
        lane += 1

# Now handle SRA-style groups: _2 -> R1, _3 -> R2, optional _1 -> I1
for group in sra_groups:
    i1 = group["i1"]
    r1 = group["r1"]
    r2 = group["r2"]

    if i1 is not None:
        i1.rename(fastq_all / f"{sample_id}_S1_L{lane:03d}_I1_001.fastq.gz")

    r1.rename(fastq_all / f"{sample_id}_S1_L{lane:03d}_R1_001.fastq.gz")
    r2.rename(fastq_all / f"{sample_id}_S1_L{lane:03d}_R2_001.fastq.gz")
    lane += 1

if lane == 1:
    # Neither bcl2fastq-style R1/R2 nor SRA-style *_1/_2/_3 could be recognised
    raise AssertionError(
        dedent(
            "Could not recognise FASTQ naming pattern as either bcl2fastq R1/R2 "
            "or SRA *_1/_2/_3. Please check filenames or open an issue on GitHub."
        )
    )

# fmt: off
run(
    [
        "cellranger", "count",
        "--id", "${prefix}",
        "--fastqs", str(fastq_all),
        "--transcriptome", "${reference.name}",
        "--localcores", "${task.cpus}",
        "--localmem", "${task.memory.toGiga()}",
        *shlex.split("""${args}"""),
    ],
    check=True,
)
# fmt: on

# Output version information
version = run(
    ["cellranger", "-V"],
    text=True,
    check=True,
    capture_output=True,
).stdout.replace("cellranger cellranger-", "")

# alas, no `pyyaml` pre-installed in the cellranger container
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'  cellranger: "{version}"\\n')
