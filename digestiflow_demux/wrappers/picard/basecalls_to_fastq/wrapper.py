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

pushd {snakemake.params.output_dir}

for sheet in {snakemake.input.sheets}; do
    prep_dir=$(dirname $sheet)
    lane=$(basename $prep_dir)

    echo $lane
    echo $sheet

    head -n 1000 $sheet

    picard IlluminaBasecallsToFastq \
        BASECALLS_DIR={snakemake.params.input_dir}/Data/Intensities/BaseCalls \
        READ_STRUCTURE={snakemake.params.read_structure} \
        LANE=$lane \
        MULTIPLEX_PARAMS=$sheet \
        BARCODES_DIR=$prep_dir \
        RUN_BARCODE={snakemake.params.run_number} \
        FLOWCELL_BARCODE={snakemake.params.flowcell_token} \
        MACHINE_NAME={snakemake.params.machine_name} \
        NUM_PROCESSORS={snakemake.threads} \
        COMPRESS_OUTPUTS=true

done

popd
"""
)
