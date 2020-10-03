# Bioinformatic-tool_for_sam-vsf_files_manipulation
Bioinformatic tool for filtering and manipulating sam/bam alignment, bcf/vcf variant calling and fastq files.

## Python modules requirements
You need to have Python version >=3.0 and the following modules installed:
<br/>numpy
<br/>pysam

## Usage instructions

## Full command-line options (--help)
```
usage: bioinformatic_formats_tool.py [-h] -o <path-to-output> [-b] [-r REF_MATCH] [-mq <1-43>] [-rq <1-40>] [-p <0-100>] [-l READ_LEN] [-s] [-m] [-so] [-ix] [-fs] [-bf]
                                     [--fastq-filter <min-length min-quality> <min-length min-quality>] [-fc] [-ff] [-v <ref-name read-depth map-qual> [<ref-name read-depth map-qual> ...]]
                                     <input file> [<input file> ...]

Bioinformatic software tool for filtering and manipulating sam/bam alignment, bcf/vcf variant calling and fastq files; version 0.2. Author: Vadim (Dani) Dubinsky (dani.dubinsky@gmail.com)

positional arguments:
  <input file>          input file (sam/bam/vcf/bcf/fastq) or files (paired fastq)

optional arguments:
  -h, --help            show this help message and exit
  -o <path-to-output>, --output <path-to-output>
                        path to output directory
  -b, --bam-filter      apply bam filter, select from the next 5 parameters (default: retain all reads)
  -r REF_MATCH, --ref-match REF_MATCH
                        bam filter, reads mapped to specific reference, wildcards (*) and regexp can be used
  -mq <1-43>, --map-qual <1-43>
                        bam filter, reads above a mapping quality (MAPQ) threshold
  -rq <1-40>, --read-qual <1-40>
                        bam filter, reads above a mean Phred score threshold
  -p <0-100>, --pident <0-100>
                        bam filter, reads above percent identity between the reads and reference
  -l READ_LEN, --read-len READ_LEN
                        bam filter, reads above a length threshold
  -s, --mapread-stats   Calculate reads mapping per each reference contig/gene
  -m, --mapped-reads    extract only mapped reads (for sam/bam files)
  -so, --sort           samtools sort command (for sam/bam files)
  -ix, --index          samtools index command (for sam/bam files)
  -fs, --flagstat       samtools flagstat command (for sam/bam files, output reads mapping summary stats)
  -bf, --bamtofastq     convert a bam to fastq file (paired reads interleaved in a single file)
  --fastq-filter <min-length min-quality> <min-length min-quality>
                        filter a single fastq or two paired fastq file/s based on read length and mean Phred score quality
  -fc, --fastq-concat   convert 2 paired-end fastq files to a single interleaved file
  -ff, --fastq-fasta    convert fastq file to fasta format
  -v <ref-name read-depth map-qual> [<ref-name read-depth map-qual> ...], --vcf-filter <ref-name read-depth map-qual> [<ref-name read-depth map-qual> ...]
                        Filter vcf/bcf variant calling file based on reference name, read depth and mapping quality
```

## Examples
