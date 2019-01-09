"""Snakemake wrapper for MultiQC.

This file is parth of Digestifly Demux.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x
set -euo pipefail

rm -f {snakemake.output}

multiqc \
    --zip-data-dir \
    --outdir $(dirname {snakemake.output.html}) \
    --interactive \
    {snakemake.input}

"""
)
