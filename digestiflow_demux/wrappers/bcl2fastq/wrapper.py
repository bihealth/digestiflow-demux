"""Snakemake wrapper for bcl2fastq2.

This file is part of Digestify Demux.
"""

import os
import sys

from snakemake import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from digestiflow_demux.bases_mask import return_bases_mask  # noqa: E402

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# Get bases mask for the current sample sheet
planned_reads = snakemake.config["flowcell"]["planned_reads"]  # noqa: F821
demux_reads = snakemake.config["flowcell"].get("demux_reads")  # noqa: F821
bases_mask_illumina = return_bases_mask(planned_reads, demux_reads)

shell.executable("/bin/bash")

# Get number of barcode mismatches, defaults to 0 for bcl2fastq v1.
barcode_mismatches = snakemake.config.get("barcode_mismatches")  # noqa
if barcode_mismatches is None:
    barcode_mismatches = 0

# More than 8 threads will not work for bcl2fastq.
bcl2fastq_threads = min(8, snakemake.config["cores"])  # noqa


shell(
    r"""
set -euo pipefail
set -x

# -------------------------------------------------------------------------------------------------
# Setup Auto-cleaned Temporary Directory.

export TMPDIR=$(mktemp -d)
mkdir -p $TMPDIR

# -------------------------------------------------------------------------------------------------
# Print Sample Sheet

head -n 10000 {snakemake.input.sheet}

# -------------------------------------------------------------------------------------------------
# Fixup Perl & Conda

# NB: there is no good Perl ecosystem on conda and dependencies are broken. We need to fix the
# installation of XML::Simple.
inc=$(perl -e "print qq(@INC)" | cut -d ' ' -f 4)
if [[ ! -e $inc/XML ]]; then
    ln -sr $(find $(dirname $inc) -name XML | grep -v x86 | sort | tail -n 1) $inc/XML
fi

# -------------------------------------------------------------------------------------------------
# Run blc2fastq v1

# Prepare output directory with Makefile etc.

configureBclToFastq.pl \
    --mismatches {barcode_mismatches} \
    --sample-sheet {snakemake.input.sheet} \
    --input-dir {snakemake.params.input_dir}/Data/Intensities/BaseCalls \
    --output-dir $TMPDIR/demux_out \
    --fastq-cluster-count 0 \
    --use-bases-mask {bases_mask_illumina} \
    --force \
    {snakemake.params.tiles_arg}

# Actually perform the demultiplexing using Make
make \
    -C $TMPDIR/demux_out \
    -j {bcl2fastq_threads}

# -------------------------------------------------------------------------------------------------
# Move Files to Destination.

# Move sample FASTQ files.
flowcell={snakemake.params.flowcell_token}
srcdir=$TMPDIR/demux_out/Project_Project

for path in $srcdir/*/*.fastq.gz; do
    sample=$(basename $(dirname $path) | cut -d _ -f 2-)
    lane=$(basename $path | rev | cut -d _ -f 3 | rev)
    dest={snakemake.params.output_dir}/$sample/$flowcell/$lane/$(basename $path)

    cp $path $dest
    pushd $(dirname $dest)
    md5sum $(basename $dest) >$(basename $dest).md5
    popd
done

# Move undetermined FASTQ files.
srcdir=$TMPDIR/demux_out

for path in $srcdir/Undetermined_indices/*/*.fastq.gz; do
    lane=$(basename $path | rev | cut -d _ -f 3 | rev)
    dest={snakemake.params.output_dir}/Undetermined/$flowcell/$lane/$(basename $path)

    cp $path $dest
    pushd $(dirname $dest)
    md5sum $(basename $dest) >$(basename $dest).md5
    popd
done

touch {snakemake.output.marker}
"""
)
