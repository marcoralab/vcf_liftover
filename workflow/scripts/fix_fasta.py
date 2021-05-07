#!/usr/bin/env python

import gzip
import re
import pandas as pd
import tempfile
import os
import shutil
import sys
from build_helpers import *


def process_line(line, build, regex, b37_map):
    line = line.strip()
    hg38 = ['GRCH38', 'B38', 'HG38', 'HG19']
    if line[0] == '>' and build.upper() not in hg38:
        contig = re.search(regex, line).group()
        if build.upper() == 'B37':
            contig = b37_convert(contig, b37_map, build.upper())
        else:
            contig = generic_convert(contig)
        return re.sub(regex, contig, line)
    return line

if __name__ == "__main__":
    if snakemake:
        fasta_file = snakemake.input['fasta']
        output_file = snakemake.output[0]
        b37_file = snakemake.input['b37']
        build = snakemake.wildcards['tobuild']
        smk = True
    else:
        fasta_file = 'data/ref/hg19.fa.gz'
        output_file = 'data/ref/b37.fa.gz'
        b37_file = 'data/ref/b37.builds.tsv'
        build = 'b37'
        smk = False

    hg38 = ['GRCH38', 'B38', 'HG38', 'HG19']

    if build in hg38:
        shutil.copy(fasta_file, output_file)
        res = os.system('samtools faidx {ofi}'.format(ofi=output_file))
        assert (res == 0), "Failed to faidx"
        sys.exit(0)

    b37 = pd.read_csv(b37_file, sep='\t')

    regex = re.compile(r'(?<=^>)\S+')

    with gzip.open(fasta_file, 'rt') as infile:
        with tempfile.NamedTemporaryFile('wt') as tfile:
            for line in infile:
                oline = process_line(line, build, regex, b37)
                print(oline, file=tfile)
            res = os.system('cat {tf} | bgzip > {ofi}'.format(
                tf=tfile.name, ofi=output_file))
            assert (res == 0), "Failed to bgzip"
            res = os.system('samtools faidx {ofi}'.format(ofi=output_file))
            assert (res == 0), "Failed to faidx"
