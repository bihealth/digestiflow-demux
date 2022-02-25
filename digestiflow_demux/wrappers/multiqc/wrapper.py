"""Snakemake wrapper for MultiQC.

This file is parth of Digestifly Demux.
"""

import tempfile

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

with tempfile.NamedTemporaryFile("wt+") as listf:
    print("\n".join(snakemake.input), file=listf)  # noqa
    listf.flush()
    listf.seek(0)

    shell(
        r"""
  set -x
  set -euo pipefail

  rm -f {snakemake.output}

  multiqc \
      --zip-data-dir \
      --outdir $(dirname {snakemake.output.html}) \
      --interactive \
      --file-list {listf.name}
  """
    )
