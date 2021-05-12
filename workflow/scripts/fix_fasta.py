#!/usr/bin/env python

from xopen import xopen
import re
import pandas as pd
import tempfile
import os
import shutil
import sys
try:
    from build_helpers import *
except:
    from workflow.scripts.build_helpers import *


def process_line(line, build, regex, b37_map):
    line = line.decode()
    contig = re.search(regex, line).group()
    b37 = ['B37', 'HUMANG1KV37', 'GRCH37']
    if build.upper() in b37:
        contig = b37_convert(contig, b37_map, build.upper())
    else:
        contig = generic_convert(contig)
    return re.sub(regex, contig, line).encode('ascii')

if __name__ == "__main__":
    if 'snakemake' in locals():
        fasta_file = snakemake.input['fasta']
        output_file = snakemake.output[0]
        b37_file = snakemake.input['b37']
        build = snakemake.wildcards['tobuild']
        smk = True
    else:
        fasta_file = 'temp/ref/hg19.fa.gz'
        output_file = 'data/ref/b37.fa.gz'
        b37_file = 'data/ref/b37.builds.tsv'
        build = 'b37'
        smk = False

    hg38 = ['GRCH38', 'B38', 'HG38', 'HG19']

    if build.upper() in hg38:
        shutil.copy(fasta_file, output_file)
        res = os.system('samtools faidx {ofi}'.format(ofi=output_file))
        assert (res == 0), "Failed to faidx"
        sys.exit(0)

    b37 = pd.read_csv(b37_file, sep='\t')

    regex = re.compile(r'(?<=^>)\S+')

    with xopen(fasta_file, 'rb') as infile:
        with tempfile.NamedTemporaryFile('wb') as tfile:
            print(tfile.name)
            for line in infile:
                if line[0] == 62: # 62 is ascii for '>'
                    line = process_line(line, build, regex, b37)
                tfile.write(line)
            res = os.system('cat {tf} | bgzip > {ofi}'.format(
                tf=tfile.name, ofi=output_file))
            assert (res == 0), "Failed to bgzip"
            res = os.system('samtools faidx {ofi}'.format(ofi=output_file))
            assert (res == 0), "Failed to faidx"
