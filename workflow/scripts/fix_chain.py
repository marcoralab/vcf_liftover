#!/usr/bin/env python

import gzip
import pandas as pd
try:
    from build_helpers import *
except:
    from workflow.scripts.build_helpers import *


def process_line(line, in_build, out_build, b37_map):
    cols = line.strip().split(' ')
    b37 = ['B37', 'HUMANG1KV37', 'GRCH37']
    if cols[0] == 'chain':
        if in_build.upper() in b37:
            cols[2] = b37_convert(cols[2], b37_map, in_build.upper())
        elif in_build.upper() not in ['GRCH38', 'B38', 'HG38', 'HG19']:
            cols[2] = generic_convert(cols[2])

        if out_build.upper() in b37:
            cols[7] = b37_convert(cols[7], b37_map, out_build.upper())
        elif out_build.upper() not in ['GRCH38', 'B38', 'HG38', 'HG19']:
            cols[7] = generic_convert(cols[7])

        return ' '.join(cols)
    else:
        return line.strip()

if __name__ == "__main__":
    if 'snakemake' in locals():
        input_file = snakemake.input['chain']
        output_file = snakemake.output[0]
        b37_file = snakemake.input['b37']
        in_build = snakemake.wildcards['frombuild']
        out_build = snakemake.wildcards['tobuild']
        smk = True
    else:
        input_file = 'temp/ref/hg38_to_hg19.over.chain.gz' #'temp/ref/hg19_to_hg38.over.chain.gz'
        output_file = 'data/ref/hg38_to_b37.over.chain.gz' #'data/ref/b37_to_b38.over.chain.gz'
        b37_file = 'data/ref/b37.builds.tsv'
        in_build = 'hg38' #'b37'
        out_build = 'b37' #'HG38'
        smk = False

    b37 = pd.read_csv(b37_file, sep='\t')
    with gzip.open(input_file, 'rt') as infile:
        with gzip.open(output_file, 'wt') as outfile:
            for line in infile:
                oline = process_line(line, in_build, out_build, b37)
                print(oline, file=outfile)
