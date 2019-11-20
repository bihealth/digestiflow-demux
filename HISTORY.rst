===========================
Digestiflow Demux Changelog
===========================

-----------------
HEAD (unreleased)
-----------------

- Adjustment to NovaSeq firmware.

------
v0.4.2
------

- Fixing interpretation of config vs. args.

------
v0.4.1
------

- Adding ``--config`` argument.
- Fixing bug with XML parsing.

------
v0.4.0
------

- Automatically selecting Picard when necessary.
- Fixing wrong assertion that disallowed RTAv3 for Picard.
- Fixing PE demultiplexing with RTAv1 (important for MiSeq).
- Fixing setting of temporary directory.
- Bumping FastQC dependency.
- Allowing for using Snakemake DRMAA execution.
- Posting Illumina demultiplexing HTML report to Digestiflow Server.
- Allowing latest RTA v1.x to be demultiplexed with bcl2fastq v2.

------
v0.3.1
------

- Support for RTAv3 (via ``bcl2fastq2``).
- Skipping libraries without a lane assignment.
- Fixing creation of BCL tarballs.
- Ignore libraries without configured lanes.
- Fix demultiplexing with bcl2fastq v1 (MiSeq).

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
