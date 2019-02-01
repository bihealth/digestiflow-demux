"""Snakemake wrapper for picard ExtractIlluminaBarcodes.

This file is part of Digestify Demux.
"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bihealth.de>"

shell.executable("/bin/bash")

# Consider tiles to process.
tiles = snakemake.config.get("tiles")  # noqa
if not tiles:
    tiles = ""
else:
    tiles = "TILE_LIMIT=1"

shell(
    r"""
set -euo pipefail
set -x

# -------------------------------------------------------------------------------------------------
# Setup Auto-cleaned Temporary Directory.

export TMPDIR=$(mktemp -d)
mkdir -p $TMPDIR

pushd $TMPDIR

for sheet in {snakemake.input.sheets}; do
    prep_dir=$(dirname $sheet)
    lane=$(basename $prep_dir)
    mkdir -p $lane

    head -n 1000 $sheet

    picard -Xmx16g \
        IlluminaBasecallsToFastq \
        BASECALLS_DIR={snakemake.params.input_dir}/Data/Intensities/BaseCalls \
        READ_STRUCTURE={snakemake.params.read_structure} \
        LANE=$lane \
        MULTIPLEX_PARAMS=$sheet \
        BARCODES_DIR=$prep_dir \
        RUN_BARCODE={snakemake.params.run_number} \
        FLOWCELL_BARCODE={snakemake.params.flowcell_token} \
        MACHINE_NAME={snakemake.params.machine_name} \
        NUM_PROCESSORS={snakemake.threads} \
        COMPRESS_OUTPUTS=true {tiles}

    # Move files to destination
    for path in $lane/*.fastq.gz; do
        sample=$(basename ${{path%%.fastq.gz}})
        read=${{sample##*.}} # 1, 2, barcode_1, barcode_2
        sample=${{sample%%.*}}
        flowcell={snakemake.params.flowcell_token}
        lane_name=L$(printf "%03d" $lane)

        # Add R to read files 1 and 2
        if [[ ${{#read}} -eq 1 ]]; then
            read=R${{read}}
        fi

        dest={snakemake.params.output_dir}/$sample/$flowcell/$lane_name/${{sample}}_${{lane_name}}_${{read}}_001.fastq.gz
        cp -dR $path $dest
        pushd $(dirname $dest)
        md5sum $(basename $dest) >$(basename $dest).md5
        popd

    done
done

popd
"""
)
