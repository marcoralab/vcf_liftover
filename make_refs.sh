#!/usr/bin/env bash
snakemake --profile lsf resources/ref/b36_to_b37.over.chain.gz resources/ref/b36_to_b38.over.chain.gz resources/ref/b37.fa.gz resources/ref/b37_to_b38.over.chain.gz resources/ref/b38.fa.gz resources/ref/b38_to_b37.over.chain.gz resources/ref/hg38_to_b37.over.chain.gz
