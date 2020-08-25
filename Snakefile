from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

startpath = '/sc/arion/projects/LOAD/Public_Data/HGDP'
startfile = '/hgdp_wgs.20190516.full.chr{chrom}.vcf.gz'

localrules: all, download_chain, download_fasta

rule all:
  input:
    expand('hgdp_wgs.20190516.PASS.GRCh{build}.chr{chrom}.vcf.gz',
           chrom = range(1,23), build = [37, 38])

rule filter:
  input: startpath + startfile
  output:
    vcf = 'hgdp_wgs.20190516.PASS.GRCh38.chr{chrom}.vcf.gz',
    tbi = 'hgdp_wgs.20190516.PASS.GRCh38.chr{chrom}.vcf.gz.tbi'
  conda: 'hgdpenv.yaml'
  shell: '''
bcftools view -f PASS {input} -Oz -o {output.vcf}
bcftools index -t {output.vcf}
'''

liftover_dir = ('ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/')

rule download_chain:
  input:
    chain = FTP.remote(liftover_dir + 'GRCh38_to_GRCh37.chain.gz'),
    checksum = FTP.remote(liftover_dir + 'CHECKSUMS')
  output: 'GRCh38_to_GRCh37.chain.gz'
  shell: '''
cp {input.chain} {output}
md5sum -c <(cat {input.checksum} | grep 'GRCh38_to_GRCh37.chain.gz')
'''

fasta = ('https://github.com/broadinstitute/gatk/raw/master/src/test/'
         'resources/large/Homo_sapiens_assembly38.fasta.gz')

rule download_fasta:
  input: HTTP.remote(fasta, allow_redirects=True),
  output: 'Homo_sapiens_assembly38.fasta'
  shell: '''
gzip -c {input} > {output}
'''

rule liftover:
  input:
    vcf = rules.filter.output.vcf,
    chain = rules.download_chain.output,
    fasta = rules.download_fasta.output
  output:
    vcf = 'hgdp_wgs.20190516.PASS.GRCh37.chr{chrom}.vcf.gz',
    tbi = 'hgdp_wgs.20190516.PASS.GRCh37.chr{chrom}.vcf.gz.tbi'
  conda: 'hgdpenv.yaml'
  shell: '''
CrossMap vcf {input.chain} {input.vcf} {input.fasta} {output.vcf}
bcftools index -t {output.vcf}
'''
