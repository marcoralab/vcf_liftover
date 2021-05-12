from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from urllib.request import urlopen
from urllib.error import URLError

try:
  response = urlopen('https://www.google.com/', timeout=10)
  iconnect = True
except urllib.error.URLError as ex:
  iconnect = False

class dummyprovider:
  def remote(string_, allow_redirects = "foo"):
    return string_

FTP = FTPRemoteProvider() if iconnect else dummyprovider
HTTP = HTTPRemoteProvider() if iconnect else dummyprovider


localrules: download_chain, download_fasta
gatk = 'docker://broadinstitute/gatk:4.1.8.1'

#liftover_dir = ('ftp://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/')


rule download_chain:
  input:
    chain = HTTP.remote('http://hgdownload.cse.ucsc.edu/goldenpath/hg{frombuild_raw_chain}/liftOver/hg{frombuild_raw_chain}ToHg{tobuild_raw_chain}.over.chain.gz', allow_redirects=True),
    checksum = HTTP.remote('http://hgdownload.cse.ucsc.edu/goldenpath/hg{frombuild_raw_chain}/liftOver/md5sum.txt', allow_redirects=True)
  params:
    chainfile = 'hg{frombuild_raw_chain}ToHg{tobuild_raw_chain}.over.chain.gz'
  output: temp('temp/ref/hg{frombuild_raw_chain}_to_hg{tobuild_raw_chain}.over.chain.gz')
  shell: '''
md5sum -c <(cat {input.checksum} | \
            grep '{params.chainfile}' | \
            sed 's|{params.chainfile}|{input.chain}|g')
cp {input.chain} {output}
'''


def calculate_build_chain_(in_build):
  if in_build.lower() in ['grch35', 'b35', 'ncbi35', 'hg17']:
    build = '17'
  if in_build.lower() in ['grch36', 'b36', 'ncbi36', 'hg18']:
    return '18'
  elif in_build.lower() in ['grch37', 'b37', 'hg19', 'humang1kv37']:
    return '19'
  elif in_build.lower() in ['grch38', 'b38', 'hg38']:
    return '38'

def calculate_build_chain(wildcards):
  fbr = calculate_build_chain_(wildcards.frombuild)
  tbr = calculate_build_chain_(wildcards.tobuild)
  return 'temp/ref/hg{bfrom}_to_hg{bto}.over.chain.gz'.format(
    bfrom=fbr, bto=tbr)

rule fix_chain:
  input:
    chain = calculate_build_chain,
    b37 = 'data/ref/b37.builds.tsv'
  output: 'data/ref/{frombuild}_to_{tobuild}.over.chain.gz'
  conda: 'envs/hgdpenv.yaml'
  script: 'scripts/fix_chain.py'


rule download_fasta:
  input:
    fasta = HTTP.remote('https://hgdownload.soe.ucsc.edu/goldenPath/{tobuild_raw}/bigZips/{tobuild_raw}.fa.gz', allow_redirects=True),
    checksum = HTTP.remote('https://hgdownload.soe.ucsc.edu/goldenPath/{tobuild_raw}/bigZips/md5sum.txt', allow_redirects=True)
  params:
    fastafile = '{tobuild_raw}.fa.gz'
  output: temp('temp/ref/{tobuild_raw,hg19|hg38}.fa.gz')
  conda: 'envs/hgdpenv.yaml'
  shell: '''
md5sum -c <(cat {input.checksum} | \
            grep '{params.fastafile}' | \
            sed 's|{params.fastafile}|{input.fasta}|g')
zcat {input.fasta} | bgzip > {output}
'''

def calculate_build_fasta(wildcards):
  build = wildcards.tobuild
  if build.lower() in ['grch35', 'b35', 'ncbi35', 'hg17']:
    build = 'hg17'
  if build.lower() in ['grch36', 'b36', 'ncbi36', 'hg18']:
    build = 'hg18'
  elif build.lower() in ['grch37', 'b37', 'hg19', 'humang1kv37']:
    build = 'hg19'
  elif build.lower() in ['grch38', 'b38', 'hg38']:
    build = 'hg38'
  return 'temp/ref/{}.fa.gz'.format(build)

rule fix_fasta:
  input:
    fasta = calculate_build_fasta,
    b37 = 'data/ref/b37.builds.tsv'
  output: 'data/ref/{tobuild}.fa.gz'
  conda: 'envs/hgdpenv.yaml'
  script: 'scripts/fix_fasta.py'

rule dict_fasta:
  input: rules.fix_fasta.output
  output: 'data/ref/{tobuild}.dict'
  container: gatk
  shell: 'gatk CreateSequenceDictionary -R {input}'
