===========================
Digestiflow Demux Changelog
===========================

-----------------
HEAD (unreleased)
-----------------

- Enable use of bcl2fastq packages from conda.

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
