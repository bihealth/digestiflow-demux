===========================
Digestiflow Demux Changelog
===========================

-----------------
HEAD (unreleased)
-----------------

- Support for RTAv3 (via ``bcl2fastq2``).
- Skipping libraries without a lane assignment.
- Fixing creation of BCL tarballs.

------
v0.3.0
------

- Enable use of bcl2fastq packages from conda.
- Parse and check custom bases masks.
- Allow passing bcl2fastq2 parameters from digestiflow-server, e.g. for 10x demultiplexing.
- Add index reads to output files.
- Remove old sample sheets before running to ensure consistency.

------
v0.2.0
------

- Also writing out log files into target directories.
- Adding ``--force-demultiplexing`` and ``--only-post-message`` command line arguments.
- Wrapped Snakemake call writes ``.gz``-compressed logs now.
- Fixing problem with synchronous reading/writing of stdout/stderr.
- Adding support for ``demux_reads`` value from Digestiflow Web.
- Adding support for demultiplexing with Picard.

----
v0.1
----

- Everything is new.
