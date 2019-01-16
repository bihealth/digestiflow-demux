|Bioconda|
|Build Status|

.. |Bioconda| image:: https://img.shields.io/conda/dn/bioconda/digestiflow-demux.svg?label=Bioconda
   :target: https://bioconda.github.io/recipes/digestiflow-demux/README.html
.. |Build Status| image:: https://travis-ci.org/bihealth/digestiflow-demux.svg?branch=master
   :target: https://travis-ci.org/bihealth/digestiflow-demux

=================
DigestiFlow Demux
=================

A command line client tool to perform semiautomatic demultiplexing of Illumina flowcells using data from Digestiflow Web.

------------
Installation
------------

The recommended way to install is to install the `digestiflow-demux` package from Bioconda.

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

    $ digestiflow-demux \
        --project-uuid UUID \
        path/to/output/dir \
        path/to/input/runs \
        path/to/input/second_run

The call will create one output directories in ``path/to/output/dir`` for each input run.
Flow cell sample sheets will be pulled from the project with the given ``UUID``.

-------------
Configuration
-------------

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
