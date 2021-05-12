# A GATK Picard-based pipeline for automatically lifting over vcf files

This pipeline does the following:

1.  Download and process reference and chain files from UCSC
2.  Split input files into chromosomes if it is whole-genome
3.  Run GATK Picard liftover
4.  Separate out variants that moved to a different chromosome
5.  Concatenate the chromosomes if requested.

These features are not yet implemented but can be added on request:

1.  Support for PLINK or oxford files
2.  Automatic contig detection
3.  Filtering of VCFs beyond using the filter column

The pipeline will parallelize where possible and can process more than one fileset simultaniously.

## Configuration file

The configuration file is config/config.yaml here are the options:

### `output_builds`:

A list of target builds for liftover. "hg", "GRCh" and "NCBI" builds will be converted to "b" builds.

### `all_GATK_builds`:

Set this to `true` or `false`.

Many b37 files are mislabled as hg19 or GRCh37, which causea liftover to fail. This option prevents that by assuming every hg19/b37/GRCh37/HumanG1Kv37 file is b37 or HumanG1Kv37.

If you are sure you have the right build, you can set this to false, otherwise, you can generally assume that files labeled as hg19 or GRCh37 are actually HumanG1Kv37 or b37, which are identical other than the inclusion of an hsv dummy contig (NC_007605). If chromosme 1 is "1", it is one of those, and if chromosme 1 is "chr1", it is hg19 or GRCh37.

GATK have an article on the discrepencies between these builds [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies).

### `concatenate`:

Set this to `true` or `false`. Concatenate the chromosomes after lifover is done.

### `inputs`:

A list of input file\[sets\] and options in the following format:

```yaml
inputs:
  output_name:
    build: build_name
    input_vcf: vcf_file_name
    contigs: 'from:to,other'
    filter: vcf_filter_string
```

#### `output_name`:

Replace `output_name` with the output stem for your file. The pipeline will create the following files:

*   `logs/{output_name}.from-{source build)_to-{destination build}.chr{chrom}.liftover.log`: The GATK liftover log file.
*   `output/{output_name}.{destination build}.chr{chrom}_only.vcf.gz`: Per-chromosome lifted-over vcf outputs where the cromosome matches the original chromosome before liftover.
*   `output/{output_name}.mismatched_chroms.vcf.gz`: A concatenated lifted-over vcf output the cromosome has changed from the original chromosome before liftover.
*   `output/{output_name}.{destination build}.same_chr.vcf.gz`: Concatenated lifted-over vcf outputs where the cromosome matches the original chromosome before liftover. Only produced if `concatenate` is set to `true`.

#### `build`:

Replace `build_name` with the original build of your file. Keep in mind that the chromosome contigs for any build starting with "hg" or "GRCh" will look like "chrN" rather than "N". We have relaxed this for hg18 and expect "N", but expect compliant contig names for builds 37 and up. for more details, see [`all_GATK_builds`](###`all_GATK_builds`:).

Valid build names include:
*   *b35:* NCBI35, hg17, b35
*   *b36:* NCBI36, hg18, b36
*   *b37:* GRCh37, hg19, b37, HumanG1Kv37
*   *b38:* GRCh38, hg38, b38

It will also work if you enter GRCh35 or GRCh36, though those builds do not actually exist.

#### `input_vcf`:

The `vcf.gz` input file(s). If `{chrom}` is in the provided string, the pipeline will read a different vcf file for each contig/chrom. If it is absent, the pipeline will split the input file into individual contigs.

If files are split by contig and the mitochondrial dna is included, `{chrom}` should be 'MT' instead of 'M' in the file name.

#### `contigs` (optional):

Contigs in the input file. For chromosomes 1-22, use raw numbers (e.g. 1 instead of chr1), use MT for mitochondria, and use X and Y for sex chromosomes. Use the VCF contig name for other contigs. (liftover results for other contigs will be in `output/{output_name}.mismatched_chroms.vcf.gz`)

Separate ranges with ":" and listed chromosomes with ",". You can put both in the same string.

If not specified, chromsomes 1-22 will be lifted over.

#### `filter` (optional):

If specified, prefilter the input VCFs based on the FILTER column. To keep variants with the PASS filter, for instance, set this to 'PASS'.

## Running

This pipeline uses anaconda environments and singularity containers. To run, use the following Snakemake command with any other desired command line options:

```bash
snakemake -j[number of jobs] --use-conda --use-singularity
```

If you use `snakejob`, run:

```bash
snakejob -j[number of jobs] --use-conda --use-singularity
```

Make sure you have access to singularity before running, and that you can pull docker containers with singularity.

## TODO

*   Get the pipeline to accept b38 alternate contigs (currently b37 only)
