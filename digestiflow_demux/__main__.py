# -*- coding: utf-8 -*-
"""Main entry point of the ``digestiflow-demux`` application."""

import os
import argparse
import logging
import shutil

import attr
import coloredlogs
import tempfile
import toml

from .exceptions import ApiProblemException
from .workflow import perform_demultiplexing

from ._version import get_versions


__version__ = get_versions()["version"]
del get_versions


@attr.s(frozen=True)
class DemuxConfig:
    """Configuration of demultiplexing app"""

    #: URL to Digestiflow Web API
    api_url = attr.ib()
    #: Token for Digestiflow Web API
    api_token = attr.ib(repr=False)
    #: Whether or not to write back to web API
    api_read_only = attr.ib()
    #: Force demultiplexing even if not marked as ready.
    force_demultiplexing = attr.ib()
    #: Whether to only post message for previous run.
    only_post_message = attr.ib()
    #: The project UUID to register the flowcell in
    project_uuid = attr.ib()

    #: Whether or not to keep the work directory
    keep_work_dir = attr.ib()
    #: Working directory
    work_dir = attr.ib()

    #: Select demultiplexing tool
    demux_tool = attr.ib()
    #: Select tiles
    tiles = attr.ib()
    #: Select lanes
    lanes = attr.ib()

    #: bcl2fastq2 parameters
    with_failed_reads = attr.ib()

    #: Degree of parallelism to use
    cores = attr.ib()
    #: Increase verbosity
    verbose = attr.ib()
    #: Decrease verbosity
    quiet = attr.ib()
    #: Whether or not to log the API token
    log_api_token = attr.ib()

    #: Path to temporary logging file
    log_path = attr.ib()

    @classmethod
    def build(cls, config, log_path=None):
        """Construct a new ``DemuxConfig`` object from configuration ``dict``."""
        return cls(
            api_url=config["web"]["url"],
            api_token=config["web"]["token"],
            api_read_only=config["web"]["api_read_only"],
            force_demultiplexing=config["web"]["force_demultiplexing"],
            only_post_message=config["web"]["only_post_message"],
            demux_tool=config["demux"]["demux_tool"],
            with_failed_reads=config["with_failed_reads"],
            project_uuid=config["demux"]["project_uuid"],
            keep_work_dir=config["demux"]["keep_work_dir"],
            work_dir=config["demux"]["work_dir"],
            tiles=config["demux"]["tiles"],
            lanes=config["demux"]["lanes"],
            cores=config["threads"],
            verbose=config["verbose"],
            quiet=config["quiet"],
            log_api_token=config["log_api_token"],
            log_path=log_path,
        )


def setup_logging(demux_config):
    """Setup logging based on ``demux_config``."""
    logger = logging.getLogger()
    # Clear logging handlers.
    logger.handlers = []

    # Setup logging to stderr
    if demux_config.quiet:
        coloredlogs.install(level=logging.WARN, logger=logger)
    elif demux_config.verbose:
        coloredlogs.install(level=logging.DEBUG, logger=logger)
    else:
        coloredlogs.install(level=logging.INFO, logger=logger)

    # Setup logging to temporary file for later posting
    handler = logging.FileHandler(os.path.join(demux_config.log_path, "digestiflow-demux.log"))
    formatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return handler  # so we can later flush()


def load_config():
    """Load configuration file"""
    config = {}
    for path in ("~/.digestiflowrc.toml", "~/digestiflowrc.toml"):
        if os.path.exists(os.path.expanduser(path)):
            logging.info("Loading configuration file %s, not looking further afterwards", path)
            with open(os.path.expanduser(path), "rt") as tomlf:
                config = toml.load(tomlf)
            break
        else:
            logging.debug("Configuration %s does not exist", path)
    else:
        logging.info("No configuration file found; not loading any")
    return config


def merge_config_args(config, args):
    """Merge args into configuration."""
    config.setdefault("web", {}).setdefault("url", None)
    if args.api_url:
        config.setdefault("web", {})["url"] = args.api_url
    config.setdefault("web", {}).setdefault("token", None)
    if args.api_token:
        config.setdefault("web", {})["token"] = args.api_url
    config.setdefault("web", {}).setdefault("force_demultiplexing", None)
    if args.force_demultiplexing:
        config["web"]["force_demultiplexing"] = args.force_demultiplexing
    config.setdefault("web", {}).setdefault("api_read_only", None)
    if args.api_read_only:
        config["web"]["api_read_only"] = args.api_read_only
    config.setdefault("web", {}).setdefault("only_post_message", None)
    if args.only_post_message:
        config["web"]["only_post_message"] = args.only_post_message
    config.setdefault("demux", {}).setdefault("project_uuid", None)
    if args.demux_tool:
        config["demux"]["demux_tool"] = args.demux_tool
    if args.project_uuid:
        config["demux"]["project_uuid"] = args.project_uuid
    config.setdefault("demux", {}).setdefault("keep_work_dir", False)
    if args.keep_work_dir:
        config["demux"]["keep_work_dir"] = True
    config.setdefault("threads", 1)
    if args.cores:
        config["threads"] = args.cores
    config.setdefault("verbose", False)
    if args.verbose is True:
        config["verbose"] = True
    config.setdefault("quiet", False)
    if args.quiet is True:
        config["quiet"] = True
    config.setdefault("with_failed_reads", False)
    if args.with_failed_reads is True:
        config["with_failed_reads"] = True

    config["log_api_token"] = args.log_api_token
    config.setdefault("demux", {})["tiles"] = args.tiles
    config.setdefault("demux", {})["lanes"] = args.lanes
    config.setdefault("demux", {})["work_dir"] = args.work_dir

    return config


def run(config, output_dir, input_dirs, log_handler):
    """Main entry point (after parsing command line options)."""
    marker_file = os.path.join(output_dir, "DIGESTIFLOW_DEMUX_DONE.txt")
    if os.path.exists(marker_file):
        logging.error("The output marker file %s already exists.", marker_file)
        logging.error("I'm refusing to overwrite it!")
        return 1

    logging.info("Starting digestiflow-demux")
    logging.info("Configuration is: %s", config)
    if config.log_api_token:
        logging.debug("  API token is %s", config.api_token)

    any_failure = False
    base_output_dir = os.path.realpath(output_dir)
    for input_dir in input_dirs:
        input_dir = os.path.realpath(input_dir)
        output_dir = os.path.join(base_output_dir, os.path.basename(input_dir))
        try:
            success, message, flowcell, client = perform_demultiplexing(
                config, input_dir, output_dir
            )
            if not success:
                any_failure = True
                logging.info(
                    "Demultiplexing failed or was not performed (flow cell not "
                    "registered or not ready?)"
                )
            # Write out log file.
            if message and flowcell and client:
                # Append log file to message in Digestiflow Web
                log_handler.flush()
                shutil.copy(log_handler.stream.name, os.path.join(output_dir, "log"))
                if not config.api_read_only:
                    client.message_attach(
                        flowcell["sodar_uuid"], message["sodar_uuid"], log_handler.stream
                    )
                # Truncate file after this flow cell to prevent confusion.
                logging.debug("Starting new log file for new flow cell")
                log_handler.stream.truncate()
        except ApiProblemException as e:
            logging.warning("There was an API problem for the flow cell. %s", e)
            logging.warning("Will continue but program will have non-zero return code.")
            any_failure = True

    return any_failure


def main(argv=None):
    """Main entry point (before parsing command line options).

    Will also load configuration TOML file at ``~/.digestiflowrc.toml`` (if any) and merge this
    with the settings from the arguments.
    """
    # Setup argument parser and parse command line arguments.
    parser = argparse.ArgumentParser(description="Run demultiplexing for Digestiflow")

    parser.add_argument(
        "--demux-tool",
        choices=["bcl2fastq", "picard"],
        default="bcl2fastq",
        help="Demultiplexing tool to use. Choices are Illumina's bcl2fastq(2) and Picard",
    )
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))
    parser.add_argument("--api-url", help="URL to Digestiflow Web API")
    parser.add_argument(
        "--api-token", default=False, help="API token to use for Digestiflow Web API"
    )
    parser.add_argument(
        "--log-api-token",
        action="store_true",
        help=(
            "Create log entry with API token (debug level; use only when debugging as this "
            "has security implications)"
        ),
    )
    parser.add_argument(
        "--api-read-only", action="store_true", help="Do not write/update flowcell info to database"
    )
    parser.add_argument(
        "--only-post-message", action="store_true", help="Only create the success message."
    )
    parser.add_argument(
        "--force-demultiplexing",
        action="store_true",
        help="Force demultiplexing even if flow cell not marked as ready",
    )
    parser.add_argument("--project-uuid", help="Project UUID to register flowcell for")
    parser.add_argument("--cores", type=int, help="Degree of parallelism to use")
    parser.add_argument("--verbose", action="store_true", default=None, help="Increase verbosity")
    parser.add_argument("--quiet", action="store_true", default=None, help="Decrease verbosity")
    parser.add_argument(
        "--keep-work-dir",
        default=False,
        action="store_true",
        help="Keep temporary working directory (useful only for debugging)",
    )
    parser.add_argument(
        "--work-dir", help="Specify working directory (instead of using temporary one)"
    )

    parser.add_argument("output_dir", metavar="OUT_DIR", help="Path to output directory")
    parser.add_argument(
        "input_dirs", metavar="SEQ_DIR", help="Path(s) to sequencer raw data directories", nargs="+"
    )

    group = parser.add_argument_group("Lane/Tile Selection")
    mutex = group.add_mutually_exclusive_group()
    mutex.add_argument(
        "--lane",
        type=int,
        default=[],
        action="append",
        dest="lanes",
        help=(
            "Select individual lanes for demultiplexing; default is to "
            "use all for which the sample sheet provides information; "
            "provide multiple times for selecting multiple lanes."
        ),
    )
    mutex.add_argument(
        "--tiles",
        type=str,
        default=[],
        action="append",
        dest="tiles",
        help=(
            "Select tile regex; provide multiple times for multiple "
            "regexes with bcl2fastq. Picard will use the first tile. "
            "Conflicts with --lane"
        ),
    )
    parser.add_argument("--with-failed-reads", action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args(argv)

    # Load configuration and merge with arguments.
    config = merge_config_args(load_config(), args)

    with tempfile.TemporaryDirectory() as tmpdir:
        demux_config = DemuxConfig.build(config, tmpdir)
        # Construct ``DemuxConfig`` from merge result and check the configuration values.
        if not demux_config.api_url:
            logging.error("The URL to the API has not been specified!")
            return 1
        if not demux_config.api_token:
            logging.error("The API token has not been specified.")
            return 1
        if not demux_config.project_uuid:
            logging.error("The project UUID is missing")
            return 1
        # Setup logging and launch.
        return run(demux_config, args.output_dir, args.input_dirs, setup_logging(demux_config))
