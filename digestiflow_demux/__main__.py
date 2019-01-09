# -*- coding: utf-8 -*-
"""Main entry point of the ``digestiflow-demux`` application."""

import os
import argparse
import logging

import attr
import toml


@attr.s
class DemuxConfig:
    """Configuration of demultiplexing app"""

    #: URL to Digestiflow Web API
    api_url = attr.ib()
    #: Token for Digestiflow Web API
    api_token = attr.ib()

    #: Degree of parallelism to use
    cores = attr.ib()
    #: Increase verbosity
    verbose = attr.ib()
    #: Decrease verbosity
    quiet = attr.ib()

    @classmethod
    def build(cls, config):
        """Construct a new ``DemuxConfig`` object from configuration ``dict``."""
        return cls(
            api_url=config["web"]["url"],
            api_token=config["web"]["token"],
            cores=config["threads"],
            verbose=config["verbose"],
            quiet=config["quiet"],
        )


def run(args):
    """Main entry point (after parsing command line options)."""


def main(argv=None):
    """Main entry point (before parsing command line options).

    Will also load configuration TOML file at ``~/.digestiflowrc.toml`` (if any) and merge this
    with the settings from the arguments.
    """
    parser = argparse.ArgumentParser(description="Run demultiplexing for Digestiflow")

    parser.add_argument("--api-url", help="URL to Digestiflow Web API")
    parser.add_argument("--api-token", help="API token to use for Digestiflow Web API")
    parser.add_argument("--cores", type=int, help="Degree of parallelism to use")
    parser.add_argument("--verbose", action="store_true", default=None, help="Increase verbosity")
    parser.add_argument("--quiet", action="store_true", default=None, help="Decrease verbosity")

    parser.add_argument("output_dir", metavar="OUT_DIR", help="Path to output directory")
    parser.add_argument(
        "input_dirs", metavar="SEQ_DIR", help="Path(s) to sequencer raw data directories", nargs="+"
    )

    args = parser.parse_args(argv)

    # Load configuration file, if any.
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

    # Now override configuration keys with values from arguments.
    config.setdefault("web", {}).setdefault("url", None)
    if args.api_url:
        config.setdefault("web", {})["url"] = args.api_url
    config.setdefault("web", {}).setdefault("token", None)
    if args.api_token:
        config.setdefault("web", {})["token"] = args.api_url
    config.setdefault("threads", 1)
    if args.threads:
        config["threads"] = args.cores
    config.setdefault("verbose", False)
    if args.verbose is True:
        config["verbose"] = True
    config.setdefault("quiet", False)
    if args.quiet is True:
        config["quiet"] = True

    # Construct ``DemuxConfig`` and check arguments.
    demux_config = DemuxConfig.build(config)
    if not demux_config.api_url:
        logging.error("The URL to the API has not been specified!")
    if not demux_config.api_token:
        logging.error("The API token has not been specified.")

    return run(demux_config)
