#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/scrnaseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'STAR': ['v_star.txt', r"(\S+)"],
    'Salmon': ['v_salmon.txt', r"salmon (\S+)"],
    'Kallisto': ['v_kallisto.txt', r"kallisto, version (\S+)"],
    'BUStools': ['v_bustools.txt', r"bustools (\S+)"],        
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['nf-core/scrnaseq'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = '<span style="color:#999999;\">N/A</span>'
results['Salmon'] = '<span style="color:#999999;\">N/A</span>'
results['Kallisto'] = '<span style="color:#999999;\">N/A</span>'
results['BUStools'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'nf-core/scrnaseq-software-versions'
section_name: 'nf-core/scrnaseq Software Versions'
section_href: 'https://github.com/nf-core/scrnaseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
