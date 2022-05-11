"""Snakemake wrapper for bclconvert.

This file is part of Digestify Demux.
"""

import os
import sys
import glob
from snakemake import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from digestiflow_demux.bases_mask import return_bases_mask  # noqa: E402
from digestiflow_demux.snakemake_support import build_sample_map  # noqa: E402

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

shell.executable("/bin/bash")

# Get number of barcode mismatches, defaults to 1.
barcode_mismatches = snakemake.config.get("barcode_mismatches")  # noqa
if barcode_mismatches is None:
    barcode_mismatches = 1

threads = min(16, snakemake.config["cores"])  # noqa

# Get bases mask for the current sample sheet
planned_reads = snakemake.config["flowcell"]["planned_reads"]  # noqa
bases_mask = os.path.basename(os.path.dirname(snakemake.input.sheet))  # noqa
bases_mask_illumina = return_bases_mask(planned_reads, bases_mask)
bases_mask_illumina = bases_mask_illumina.replace(",", ";")

# Get sample map to get sample numbering correct
sample_map = build_sample_map(snakemake.config["flowcell"])  # noqa
sample_map["Undetermined"] = "S0"

# TODO
# # Get additional parameters and write to section in sheet
# add_params = ""
# for arg, value in snakemake.config["bcl2fastq2_params"].items():  # noqa
#     if value:
#         if isinstance(value, bool):
#             add_params += " --" + arg.replace("_", "-")
#         elif isinstance(value, int):
#             add_params += " --" + arg.replace("_", "-") + " " + str(value)

shell(
    r"""
set -euo pipefail
set -x

# -------------------------------------------------------------------------------------------------
# Setup Auto-cleaned Temporary Directory.

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# -------------------------------------------------------------------------------------------------
# Manipulate and print Sample Sheet

echo "[Settings]" >> {snakemake.input.sheet}
echo "OverrideCycles,{bases_mask_illumina}" >> {snakemake.input.sheet}

head -n 10000 {snakemake.input.sheet}


# -------------------------------------------------------------------------------------------------
# Execute bclconvert

bclconvert \
    --sample-sheet {snakemake.input.sheet} \
    --bcl-input-directory {snakemake.params.input_dir} \
    --output-dir $TMPDIR/demux_out \
    --bcl-num-conversion-threads {threads} \
    --bcl-sampleproject-subdirectories true

    # add for debug purposes. The first-tiles-only fails for some Novaseq flowcells.
    #--first-tile-only true
    #--tiles 2101

ls -lh $TMPDIR/demux_out/*

# -------------------------------------------------------------------------------------------------
# Move Files to Destination.

# Move sample FASTQ files.
flowcell={snakemake.params.flowcell_token}
srcdir=$TMPDIR/demux_out/

for path in $srcdir/*/*.fastq.gz; do
    sample=$(basename $path | rev | cut -d _ -f 5- | rev)
    project=$(basename $(dirname $path))
    dest={snakemake.params.output_dir}/$project/$flowcell/$sample/{bases_mask}__$(basename $path)
    mkdir -p $(dirname $dest)
    cp -dR $path $dest
done

# Move undetermined FASTQ files.
srcdir=$TMPDIR/demux_out

for path in $srcdir/Undetermined_*; do
    if [[ -f $path ]]; then
        dest={snakemake.params.output_dir}/Undetermined/$flowcell/{bases_mask}__$(basename $path)
        mkdir -p $(dirname $dest)
        cp -dR $path $dest
    fi
done

# Write out the html report files as an archive
srcdir=$TMPDIR/demux_out
pushd $srcdir
tar -czvf {snakemake.params.output_dir}/html_report_{bases_mask}.tar.gz Reports
popd

touch {snakemake.output.marker}
"""
)


# Counting samples from 0 for each bases mask breaks snakemake's assumption that
# all samples from all sample sheets are numbered incrementally. Here, we look up the correct
# sample number and replace it in the path. Easier to do in python than in bash above.
fls = glob.glob(os.path.join(snakemake.params.output_dir, "*/*/*/" + "*.fastq.gz"))  # noqa
for path in fls:
    f = os.path.basename(path)
    name = f.split("__", 1)[1]  # remove $bases_mask__
    name_elements = name.split("_")[:-3]
    oldS = name_elements[-1]
    name = "_".join(name_elements[:-1])
    newS = sample_map[name]
    newpath = path.replace("_".join([name, oldS]), "_".join([name, newS]))
    os.rename(path, newpath)
