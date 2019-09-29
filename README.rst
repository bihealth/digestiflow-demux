|Bioconda|
|Build Status|

.. |Bioconda| image:: https://img.shields.io/conda/dn/bioconda/digestiflow-demux.svg?label=Bioconda
   :target: https://bioconda.github.io/recipes/digestiflow-demux/README.html
.. |Build Status| image:: https://travis-ci.org/bihealth/digestiflow-demux.svg?branch=master
   :target: https://travis-ci.org/bihealth/digestiflow-demux

=================
DigestiFlow Demux
=================

A command line client tool to perform semi automatic demultiplexing of Illumina flowcells using data from Digestiflow Web.

------------
Installation
------------

The recommended way is to install the ``digestiflow-demux`` package from Bioconda.

-----
Usage
-----

First, create a ``~/.digestiflowrc.toml`` file to configure access to the Digestiflow Web API.

::

    [web]
    # URL to your Digestiflow instance. "$url/api" must be the API entry URL.
    url = "https://flowcells.example.org"
    # The secret token to use for the the REST API, as created through the Web UI.
    token = "secretsecretsecretsecretsecretsecretsecretsecretsecretsecretsecr"

Then, you can use the ``digestiflow-demux`` command for demultiplexing Illumina sequencing runs.

::

    $ digestiflow-demux --project-uuid UUID OUT_PATH [IN_PATH ...]

The call will create one output directories in ``OUT_PATH`` for each ``IN_PATH`` run.
Flow cell sample sheets will be pulled from the project with the given ``UUID``.

Configuration
=============

The following shows all configuration settings loaded from ``~/.digestiflowrc.toml`` with their defaults.

::

    [web]
    # URL to your Digestiflow instance. "$url/api" must be the API entry URL.
    url =
    # The secret token to use for the the REST API, as created through the Web UI.
    token =

    [demux]
    # UUID of the project to import to.
    project_uuid =
    # Whether or not to keep the working directory around (useful for debugging only).
    keep_work_dir = false
    # Whether or not to increase verbosity.
    verbose = false
    # Whether or not to decrease verbosity.
    quiet = false

You can configure cluster execution in the configuration file by the following settings in the ``[demux]`` section.

::

    [demux]
    # The number of cores to use.
    threads = 32
    # DRMAA command line to use (see snakemake's --drmaa), note the leading space.
    drmaa = " -V -cwd -S /bin/bash ..."
    # Cluster configuration JSON file.
    cluster_config = "/home/demux_user/.digestiflow.cluster_config.json"
    # Optional path to Snakemake job script to use.
    jobscript = "/home/demux_user/.digestiflow.jobscript.sh"

Calling
=======

Generally, to run the demultiplexing the flow cells at ``IN_PATH`` and ``IN_PATH2`` in the project with UUID ``UUID``, use the following command:

::

    digestiflow-demux --project-uuid UUID OUT_PATH IN_PATH IN_PATH2

The command line help is available through

:::

    digestiflow-demux --help

----------------------
Demultiplexing Process
----------------------

The demultiplexing will run as follows for each input directory.

1. The target directory is a directory below the output path with the basename of the current input directory.
2. In case that a file ``DIGESTIFLOW_DEMUX_DONE.txt`` already exists, the program halts.
3. The flow cell information is queried through the Digestiflow API.
4. If the demultiplexing state is not "is ready" then this input directory is skipped.
5. The sample sheet information is queried through the API and a sample sheet is written to the output directory.
6. A working directory is created and a Snakemake workflow is prepared that will perform the actual demultiplexing.
   The working directory will be automatically removed when the program terminates.
7. The ``snakemake`` executable is called internally (a log file is kept) which will performs the necessary steps:

    a. Call Illumina ``bcl2fastq`` in the appropriate version for the flow cell RTA version.
       Optionally, Picard can be used for demultiplexing of RTA v2 files.
       Compute MD5 sums for all read files.
    b. Run FastQC over the output files.
    c. Run MultiQC for aggregating the FastQC report.

7. The status of the flow cell is updated through the REST API to success/failed.
8. A message is created through the Digestiflow REST API indicating success or failure with the attachments

    - MultiQC HTML report,
    - MultiQC data ZIP file,
    - Snakemake call log file,
    - ``digestiflow-demux`` log file

The behaviour of the program can be changed using the following arguments.

- ``--demux-tool {bcl2fastq,picard}`` -- allow selection of demultiplexing pipeline
- ``--api-read-only`` -- only get information from API, do not write
- ``--only-post-message`` -- skip demultiplexing and only post message through API with collected information
- ``--force-demultiplexing`` -- force the demultiplexing even if the flow cell is not marked as "ready"
- ``--keep-work-dir`` -- keep the working directory (otherwise deleted on program's end)
- ``--work-dir`` -- specify explicit working directory name (implies ``--keep-work-dir``)
- ``--lane`` -- select individual lanes for demultiplexing only (``--lane 1 --lane 2`` for selecting lanes 1 and 2).
- ``--tiles`` -- select individual tiles for demultiplexing only.
  A current limitation is that you have to select tiles from all lanes or the Snakemake workflow will fail because it expects output from all lanes.
  This is mutually exclusive with ``--lane``.

The remaining arguments are self-explanatory and explain logging verbosity, and number of cores to use for the demultiplexing and QC.
