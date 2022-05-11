===========================
Digestiflow Demux Changelog
===========================

------
v0.6.0-rc.1
------

- Add bclconvert wrapper.
- Introduce new project structure (bclconvert only)

------
v0.5.4
------

- Fixing ``setup.py`` problem with 0.5.4.

------
v0.5.3
------

- Making multiqc wrapper work even for many files.
- Proper parsing for full RTA version.
- Fixing submission in case of non-created bcl2fast2 log files.
- Improving splitting heuristic.
- Allowing user to pick conda frontend.
- Fixing problem with too many output files.
- Small updates and fixes.

------
v0.5.2
------

- Adding the possibility to filter input folders by name to the flow cells that are marked as ready on the server.

------
v0.5.1
------

- Do not attempt to process empty flow cell.

------
v0.5.0
------

- Allowing to run without cluster.
- Adjustment to NovaSeq firmware.
- Fixing base mask computation.

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
