"""Snakemake wrapper for picard ExtractIlluminaBarcodes.

This file is part of Digestify Demux.
"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bihealth.de>"

shell.executable("/bin/bash")

# Get number of barcode mismatches, defaults to 0 for bcl2fastq v1.
barcode_mismatches = snakemake.config.get("barcode_mismatches")
if barcode_mismatches is None:
    barcode_mismatches = 0


shell(
    r"""
set -euo pipefail
set -x

# -------------------------------------------------------------------------------------------------
# Setup Auto-cleaned Temporary Directory.

export TMPDIR=$(mktemp -d)
mkdir -p $TMPDIR

# -------------------------------------------------------------------------------------------------

picard ExtractIlluminaBarcodes \
    BARCODE_FILE={snakemake.input} \
    BASECALLS_DIR={snakemake.params.input_dir}/Data/Intensities/BaseCalls \
    LANE={snakemake.wildcards.lane} \
    METRICS_FILE={snakemake.output} \
    OUTPUT_DIR=$(dirname {snakemake.output}) \
    READ_STRUCTURE={snakemake.params.read_structure}

"""
)
