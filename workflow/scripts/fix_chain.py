#!/usr/bin/env python

import gzip
import pandas as pd
from build_helpers import *


def process_line(line, in_build, out_build, b37_map):
    cols = line.strip().split(' ')
    if cols[0] == 'chain':
        if in_build.upper() == 'B37':
            cols[2] = b37_convert(cols[2], b37_map, in_build)
        elif in_build.upper() not in ['GRCH38', 'B38', 'HG38', 'HG19']:
            cols[2] = generic_convert(cols[2])

        if out_build.upper() == 'B37':
            cols[7] = b37_convert(cols[7], b37_map, in_build)
        elif out_build.upper() not in ['GRCH38', 'B38', 'HG38', 'HG19']:
            cols[7] = generic_convert(cols[7])

        return ' '.join(cols)
    else:
        return line.strip()

if __name__ == "__main__":
    if snakemake:
        input_file = snakemake.input['chain']
        output_file = snakemake.output[0]
        b37_file = snakemake.input['b37']
        in_build = snakemake.wildcards['frombuild']
        out_build = snakemake.wildcards['tobuild']
        smk = True
    else:
        input_file = 'data/ref/hg19_to_hg38.over.chain.gz'
        output_file = 'data/ref/b37_to_b38.over.chain.gz'
        b37_file = 'data/ref/b37.builds.tsv'
        in_build = 'b37'
        out_build ='HG38'
        smk = False

    b37 = pd.read_csv(b37_file, sep='\t')
    with gzip.open(input_file, 'rt') as infile:
        with gzip.open(output_file, 'wt') as outfile:
            for line in infile:
                oline = process_line(line, in_build, out_build, b37)
                print(oline, file=outfile)
