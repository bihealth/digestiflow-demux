"""Snakemake wrapper for FastQC.

This file is parth of Digestifly Demux.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

shell(
    r"""
fastqc \
    -o $(dirname {snakemake.output.html}) \
    {snakemake.input.fastq}

cd $(dirname {snakemake.output.html})

for x in *.html *.zip; do
    md5sum $x >$x.md5
done
"""
)
