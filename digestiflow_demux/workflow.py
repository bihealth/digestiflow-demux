"""Implementation of the workflow for demultiplexing sequencing directories."""

import csv
import json
import logging
import os
import tempfile
import xml.etree.ElementTree as ET

import snakemake
from snakemake.exceptions import WorkflowError

from digestiflow_demux import __version__
from .api_client import ApiClient, ApiException
from .exceptions import ApiProblemException

#: Path to the Snakefile.
PATH_SNAKEFILE = os.path.abspath(os.path.join(os.path.dirname(__file__), "Snakefile"))

#: Template for the success message.
TPL_MSG_SUCCESS = r"""
The demultiplexing succeeded for flow cell {flowcell[vendor_id]}.

See the attached files for quality reports.

-- 
This message was auto-created by digestiflow-demux v{version}.
"""

#: Template for the failure message.
TPL_MSG_FAILURE = r"""
The attempted demultiplexing for flow cell {flowcell[vendor_id]} has failed.

To try again, clean up any output files and mark as "ready" for demultiplexing again.

-- 
This message was auto-created by digestiflow-demux v{version}.
"""


def write_sample_sheet_v1(writer, flowcell, libraries):
    """Write V1 sample sheet"""
    header = [
        "FCID",
        "Lane",
        "SampleID",
        "SampleRef",
        "Index",
        "Description",
        "Control",
        "Recipe",
        "Operator",
        "SampleProject",
    ]
    writer.writerow(header)

    recipe = "PE_indexing" if flowcell["planned_reads"].count("T") > 1 else "SE_indexing"
    for lib in libraries:
        if lib["barcode2"]:
            barcode = "".join((lib["barcode"], "-", lib["barcode2"]))
        else:
            barcode = lib["barcode"]
        for lane in sorted(lib.lanes):
            data = [
                flowcell["vendor_id"],
                lane,
                lib["name"],
                lib["reference"],
                barcode,
                "",
                "N",
                recipe,
                lib.barcode_seq(1),
                "",
                "N",
                recipe,
                flowcell["operator"],
                "Project",
            ]
            writer.writerow(list(map(str, data)))


def write_sample_sheet_v2(writer, flowcell, libraries):
    """Write V2 sample sheet"""
    # Write [Data] Section
    writer.writerow(["[Data]"])
    dual_indexing = any(library["barcode2"] for library in libraries)
    if dual_indexing:
        writer.writerow(["lane", "sample_id", "index", "index2", "sample_project"])
    else:
        writer.writerow(["lane", "sample_id", "index", "sample_project"])

    rows = []
    for lib in libraries:
        for lane in sorted(lib["lanes"]):
            row = [lane, lib["name"], lib["barcode"]]
            if dual_indexing:
                row.append(lib["barcode2"])
            row.append("Project")
            rows.append(row)
    for row in sorted(rows):
        writer.writerow(list(map(str, row)))


def reverse_complement(seq):
    """Return reverse-complemented version of ``seq``."""
    mapping = {"A": "T", "a": "t", "C": "G", "c": "g", "G": "C", "g": "c", "T": "A", "t": "a"}
    return "".join(reversed([mapping.get(i, i) for i in seq]))


def load_run_info(path_run_info_xml):
    """Load information from ``RunInfo.xml`` file."""
    with open(path_run_info_xml, "rt") as xmlf:
        xmls = xmlf.read()
    root = ET.fromstring(xmls)
    tag_run = root.find("Run")
    return {
        "run_id": tag_run.attrib["Id"],
        "instrument": tag_run.find("Instrument").text,
        "run_no": tag_run.attrib["Number"],
        "flowcell": tag_run.find("Flowcell").text,
    }


def create_sample_sheet(config, input_dir, output_dir):
    """Query the Digestiflow API for the necessary information for building the sample sheet."""
    logging.info("Perform API queries and create sample sheet")
    client = ApiClient(
        api_url=config.api_url, api_token=config.api_token, project_uuid=config.project_uuid
    )

    logging.debug("Parsing RunInfo.xml file")
    run_info = load_run_info(os.path.join(input_dir, "RunInfo.xml"))

    logging.debug("Querying API for flow cell")
    try:
        flowcell = client.flowcell_resolve(
            instrument_id=run_info["instrument"],
            run_no=run_info["run_no"],
            flowcell_id=run_info["flowcell"],
        )
    except ApiException as e:
        raise ApiProblemException("Problem querying API for flow cell") from e

    if flowcell is None:
        logging.warning("Could not resolve flow cell via API. Not proceeding.")
        return None
    if flowcell["status_conversion"] != "ready":
        logging.warning('Status is not "ready", will skip flow cell.')
        return None

    try:
        client.flowcell_update(flowcell["sodar_uuid"], status_conversion="in_progress")
    except ApiException as e:
        raise ApiProblemException('Could not update conversion status to "in_progress"') from e

    logging.debug("Querying API for sequencing machine information")
    try:
        sequencer = client.sequencer_retrieve(sequencer=run_info["instrument"])
    except ApiException as e:
        raise ApiProblemException("Problem querying API for sequencer") from e

    logging.debug("Querying for barcode information")
    libraries = []
    for library in flowcell["libraries"]:
        if library.get("barcode_seq"):
            barcode_seq = library.get("barcode_seq")
        elif library.get("barcode"):
            try:
                barcode = client.barcodesetentry_retrieve(barcodesetentry=library.get("barcode"))
            except ApiException as e:
                raise ApiProblemException("Problem querying API for barcode #1") from e
            barcode_seq = barcode["sequence"]
        else:
            barcode_seq = ""
        if library.get("barcode_seq2"):
            barcode_seq2 = library.get("barcode_seq2")
        elif library.get("barcode2"):
            try:
                barcode2 = client.barcodesetentry_retrieve(barcodesetentry=library.get("barcode2"))
            except ApiException as e:
                raise ApiProblemException("Problem querying API for barcode #2") from e
            barcode_seq2 = barcode2["sequence"]
        else:
            barcode_seq2 = ""
        if sequencer["dual_index_workflow"] == "B":
            barcode_seq2 = reverse_complement(barcode_seq2)
        libraries.append(
            {
                "name": library["name"],
                "reference": library["reference"],
                "barcode": barcode_seq,
                "barcode2": barcode_seq2,
                "lanes": library["lane_numbers"],
            }
        )

    logging.debug("Writing out demultiplexing configuration")
    # Get barcode mismatch count or default.
    if flowcell["barcode_mismatches"] is None:
        if flowcell["rta_version"] == 1:
            barcode_mismatches = 0
        else:
            barcode_mismatches = 1
    else:
        barcode_mismatches = flowcell["barcode_mismatches"]
    with open(os.path.join(output_dir, "demux_config.json"), "wt") as jsonf:
        config_json = {
            "cores": config.cores,
            "rta_version": flowcell["rta_version"],
            "barcode_mismatches": barcode_mismatches,
            "input_dir": input_dir,
            "output_dir": output_dir,
            "flowcell": {**flowcell, "libraries": libraries},
            "tiles": config.tiles,
            "lanes": config.lanes,
        }
        json.dump(config_json, jsonf)

    logging.debug("Writing out sample sheet information")
    with open(os.path.join(output_dir, "SampleSheet.csv"), "wt") as csvf:
        funcs = {1: write_sample_sheet_v1, 2: write_sample_sheet_v2}
        funcs[flowcell["rta_version"]](csv.writer(csvf), flowcell, libraries)

    return flowcell  # Everything is fine


def send_flowcell_success_message(client, flowcell, output_dir):
    client.message_send(
        flowcell_uuid=flowcell["sodar_uuid"],
        subject="Demultiplexing failed for flow cell %s" % flowcell["vendor_id"],
        body=TPL_MSG_FAILURE.format(flowcell=flowcell, version=__version__),
        attachments=[
            os.path.join(output_dir, "multiqc/multiqc_report.html"),
            os.path.join(output_dir, "multiqc/multiqc_data.zip"),
        ],
    )


def send_flowcell_failure_message(client, flowcell):
    client.message_send(
        flowcell_uuid=flowcell["sodar_uuid"],
        subject="Demultiplexing failed for flow cell %s" % flowcell["vendor_id"],
        body=TPL_MSG_FAILURE.format(flowcell=flowcell, version=__version__),
    )


def launch_snakemake(config, flowcell, output_dir, work_dir):
    """Launch Snakemake and execute the demultiplexing"""
    logging.info("Temporary directory is %s", work_dir)
    logging.info("Start Snakemake workflow for demultiplexing")

    client = ApiClient(
        api_url=config.api_url, api_token=config.api_token, project_uuid=config.project_uuid
    )

    argv = [
        "--snakefile",
        PATH_SNAKEFILE,
        "--directory",
        work_dir,
        "--configfile",
        os.path.join(output_dir, "demux_config.json"),
        "--cores",
        config.cores,
        "--use-conda",
        "--config",
    ]
    if config.verbose:
        argv += ["--verbose", "--printshellcmds"]
    argv = list(map(str, argv))
    logging.info("Executing: snakemake %s", " ".join(argv))
    try:
        snakemake_args = {
            "snakefile": PATH_SNAKEFILE,
            "workdir": work_dir,
            "configfile": os.path.join(output_dir, "demux_config.json"),
            "cores": config.cores,
            "use_conda": True,
            "verbose": config.verbose,
            "printshellcmds": config.verbose,
        }
        # TODO: save log output of snakemake run and attach to message
        failure = not snakemake.snakemake(**snakemake_args)
    except WorkflowError as e:
        logging.warning("Running demultiplexing failed: %s", e)
        failure = True

    if not failure:
        send_flowcell_success_message(client, flowcell, output_dir)
        logging.info("Marking flowcell as complete...")
        try:
            client.flowcell_update(flowcell["sodar_uuid"], status_conversion="complete")
        except ApiException as e:
            logging.warning("Could not update conversion state to complete via API: %s", e)
        logging.info("Done running Snakemake.")
    elif flowcell:
        send_flowcell_failure_message(client, flowcell)
        logging.info("Marking flowcell as failed...")
        try:
            client.flowcell_update(flowcell["sodar_uuid"], status_conversion="failed")
        except ApiException as e:
            logging.warning("Could not update conversion state to failed via API: %s", e)
    return not failure


def perform_demultiplexing(config, input_dir, output_dir):
    """Prepare and execute the demultiplexing with the Snakemake workflow."""
    logging.info("Starting to process input directory %s", input_dir)
    logging.info("Output will go to %s", output_dir)

    logging.debug("Creating output directory %s", output_dir)
    os.makedirs(output_dir, exist_ok=True)

    flowcell = create_sample_sheet(config, input_dir, output_dir)
    if not flowcell:
        return False

    if config.work_dir:
        logging.info("Using work directory %s", config.work_dir)
        return launch_snakemake(config, flowcell, output_dir, config.work_dir)
    elif config.keep_work_dir:
        logging.info("Setup non-temporary work directory")
        return launch_snakemake(config, flowcell, output_dir, tempfile.mkdtemp("-cubi-demux"))
    else:
        logging.info("Setup temporary work directory")
        with tempfile.TemporaryDirectory("-cubi-demux") as work_dir:
            return launch_snakemake(config, flowcell, output_dir, work_dir)
