"""Implementation of the workflow for demultiplexing sequencing directories."""

import collections
import csv
import gzip
import itertools
import json
import logging
import os
import shutil
import subprocess
import sys
from threading import Thread, Lock
import tempfile
import xml.etree.ElementTree as ET

from snakemake.exceptions import WorkflowError

from digestiflow_demux import __version__
from .bases_mask import return_bases_mask
from .api_client import ApiClient, ApiException
from .exceptions import ApiProblemException, MissingOutputFile

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

    demux_reads = flowcell.get("demux_reads") or flowcell["planned_reads"]
    recipe = "PE_indexing" if demux_reads.count("T") > 1 else "SE_indexing"
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


def write_sample_sheets_v2(flowcell, libraries, output_dir):
    """Write V2 sample sheets. Write one sample sheet for each bases_mask in the config."""
    # re-shuffle dict from lib - lane - bases_mask to bases_mask - lib
    d = collections.defaultdict(dict)
    for key, lib in enumerate(libraries):
        d[lib.get("demux_reads_override", flowcell["demux_reads"])][key] = lib

    for bases_mask, libraries in d.items():
        os.makedirs(
            os.path.join(output_dir, "illumina_basesmask/{}".format(bases_mask)), exist_ok=True
        )
        with open(
            os.path.join(output_dir, "illumina_basesmask/{}/SampleSheet.csv".format(bases_mask)),
            "w",
        ) as f:
            writer = csv.writer(f, delimiter=",")
            write_sample_sheet_v2(writer, flowcell, libraries.values())


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
            barcodes = lib["barcode"].split(",")
            for barcode in barcodes:
                row = [lane, lib["name"], barcode]
                if dual_indexing:
                    row.append(lib["barcode2"])
                row.append("Project")
                rows.append(row)
    for row in sorted(rows):
        writer.writerow(list(map(str, row)))


def write_sample_sheet_picard(flowcell, libraries, output_dir):
    """Write picard sample sheets, one per lane."""
    dual_indexing = any(library["barcode2"] for library in libraries)
    if not dual_indexing:
        head_barcodes = ["barcode_sequence_1", "barcode_name", "library_name"]
        head_samplesheet = ["OUTPUT_PREFIX", "BARCODE_1"]
    else:
        head_barcodes = ["barcode_sequence_1", "barcode_sequence_2", "barcode_name", "library_name"]
        head_samplesheet = ["OUTPUT_PREFIX", "BARCODE_1", "BARCODE_2"]

    # re-shuffle dict from lib - lane - barcode to lane - lib - barcode because picard works on lanes
    d = collections.defaultdict(dict)
    for lib in libraries:
        for lane in sorted(lib["lanes"]):
            d[lane][lib["name"]] = lib

    # add Undetermined to samplesheet as picard crashes otherwise
    for lane in d:
        d[lane]["Undetermined"] = {"name": "Undetermined", "barcode": "N", "barcode2": ""}
        if dual_indexing:
            d[lane]["Undetermined"]["barcode2"] = "N"

    for lane, libraries in d.items():
        barcode_rows = []
        samples_rows = []
        for lib in libraries.values():
            output_prefix = "{lane}/{name}".format(
                name=lib["name"], flowcell=flowcell["vendor_id"], lane=lane
            )

            if dual_indexing:
                # we do not pass the barcodes names, so we use the sample name.
                barcode_row = [lib["barcode"], lib["barcode2"], lib["name"], lib["name"]]
                samples_row = [output_prefix, lib["barcode"], lib["barcode2"]]
            else:
                barcode_row = [lib["barcode"], lib["name"], lib["name"]]
                samples_row = [output_prefix, lib["barcode"]]

            # barcode file should not contain dummy for unbarcoded reads, but samplesheet must.
            if not lib["name"] == "Undetermined":
                barcode_rows.append(barcode_row)
            samples_rows.append(samples_row)

        os.makedirs(os.path.join(output_dir, "picard_barcodes/{}".format(lane)), exist_ok=True)

        with open(
            os.path.join(output_dir, "picard_barcodes/{}/barcodes.txt".format(lane)), "w"
        ) as bf, open(
            os.path.join(output_dir, "picard_barcodes/{}/samplesheet.txt".format(lane)), "w"
        ) as sf:
            barcodewriter = csv.writer(bf, delimiter="\t")
            sampleswriter = csv.writer(sf, delimiter="\t")

            barcodewriter.writerow(head_barcodes)
            sampleswriter.writerow(head_samplesheet)

            for row in sorted(barcode_rows):
                barcodewriter.writerow(list(map(str, row)))
            for row in sorted(samples_rows):
                sampleswriter.writerow(list(map(str, row)))


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


def create_sample_sheet(config, input_dir, output_dir):  # noqa: C901
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
    if flowcell["status_conversion"] != "ready" and not config.force_demultiplexing:
        logging.warning('Status is not "ready", will skip flow cell.')
        return None

    if not config.api_read_only:
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
    demux_reads_override = set()
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

        if library["demux_reads"]:
            demux_reads = library["demux_reads"]
        else:
            demux_reads = flowcell["demux_reads"] or flowcell["planned_reads"]
        demux_reads = return_bases_mask(flowcell["planned_reads"], demux_reads, "picard")
        demux_reads_override.add(demux_reads)

        libraries.append(
            {
                "name": library["name"],
                "reference": library["reference"],
                "barcode": barcode_seq,
                "barcode2": barcode_seq2,
                "lanes": library["lane_numbers"],
                "demux_reads_override": demux_reads,
            }
        )

    # Normalize bases masks, decide if paired-end, find all custom bases_masks
    planned_reads = flowcell["planned_reads"]
    demux_reads = flowcell.get("demux_reads") or planned_reads
    demux_reads = return_bases_mask(planned_reads, demux_reads, "picard")
    flowcell["demux_reads"] = demux_reads  # not used by bcl2fastq2
    flowcell["demux_reads_override"] = list(demux_reads_override)

    if config.demux_tool == "bcl2fastq" and flowcell["rta_version"] == 2:
        demux_tool = "bcl2fastq2"
    elif config.demux_tool == "bcl2fastq" and flowcell["rta_version"] == 1:
        demux_tool = "bcl2fastq1"
    else:
        demux_tool = "picard"

    bcl2fastq2_params = {
        "with_failed_reads": config.with_failed_reads,
        "create_fastq_for_index_reads": flowcell["create_fastq_for_index_reads"],
        "minimum_trimmed_read_length": flowcell["minimum_trimmed_read_length"],
        "mask_short_adapter_reads": flowcell["mask_short_adapter_reads"],
    }

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
            "demux_tool": demux_tool,
            "bcl2fastq2_params": bcl2fastq2_params,
        }
        json.dump(config_json, jsonf)

    logging.debug("Writing out sample sheet information")
    if demux_tool == "bcl2fastq1":
        with open(os.path.join(output_dir, "SampleSheet.csv"), "wt") as csvf:
            write_sample_sheet_v1(csv.writer(csvf), flowcell, libraries)
    elif demux_tool == "picard":
        write_sample_sheet_picard(flowcell, libraries, output_dir)
    else:
        write_sample_sheets_v2(flowcell, libraries, output_dir)

    return flowcell  # Everything is fine


def send_flowcell_success_message(client, flowcell, output_dir, *log_files):
    # Create renamed (and potentially compressed files
    path_in = os.path.join(output_dir, "multiqc/multiqc_%s")
    with tempfile.TemporaryDirectory() as tempdir:
        path_out = os.path.join(tempdir, "MultiQC_%%s_%s.%%s" % flowcell["vendor_id"])
        with open(path_in % "report.html", "rb") as f_in:
            with gzip.open(path_out % ("Report", "html.gz"), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        shutil.copyfile(path_in % "data.zip", path_out % ("Data", "zip"))

        # Post with renamed files.
        return client.message_send(
            flowcell_uuid=flowcell["sodar_uuid"],
            subject="Demultiplexing succeeded for flow cell %s" % flowcell["vendor_id"],
            body=TPL_MSG_SUCCESS.format(flowcell=flowcell, version=__version__),
            attachments=list(
                itertools.chain(
                    [path_out % ("Report", "html.gz"), path_out % ("Data", "zip")], log_files
                )
            ),
        )


def send_flowcell_failure_message(client, flowcell, *log_files):
    return client.message_send(
        flowcell_uuid=flowcell["sodar_uuid"],
        subject="Demultiplexing FAILED for flow cell %s" % flowcell["vendor_id"],
        body=TPL_MSG_FAILURE.format(flowcell=flowcell, version=__version__),
        attachments=log_files,
    )


def async_tee_pipe(process, input_file, out_file, out_file2, mutex):
    """Async tee-piping from input_file to two output files using the mutex."""
    logging_thread = Thread(target=tee_pipe, args=(process, input_file, out_file, out_file2, mutex))
    logging_thread.start()

    return logging_thread


def tee_pipe(process, input_file, out_file, out_stream, mutex):
    """Tee-piping from input_file to two output files using the mutex."""
    while 1:
        line = input_file.readline()
        if not line and process.poll() is not None:
            break
        else:
            with mutex:
                out_stream.write(line.decode("utf-8"))
                out_file.write(line)


def launch_snakemake(config, flowcell, output_dir, work_dir):
    """Launch Snakemake and execute the demultiplexing"""
    logging.info("Temporary directory is %s", work_dir)
    logging.info("Start Snakemake workflow for demultiplexing")

    client = ApiClient(
        api_url=config.api_url, api_token=config.api_token, project_uuid=config.project_uuid
    )

    output_log_dir = os.path.join(output_dir, "log")
    output_qc_dir = os.path.join(output_dir, "multiqc")
    if config.only_post_message:
        for path in (
            os.path.join(output_log_dir, "digestiflow-demux-snakemake.log.gz"),
            os.path.join(output_log_dir, "digestiflow-demux.log"),
            os.path.join(output_qc_dir, "multiqc_data.zip"),
            os.path.join(output_qc_dir, "multiqc_report.html"),
        ):
            if not os.path.exists(path):
                raise MissingOutputFile("Cannot post message with %s missing" % path)

    if config.only_post_message:
        logging.info("Only posting message, not running demultiplexing itself.")
        failure = False
    else:
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
            # Launch Snakemake
            proc = subprocess.Popen(
                ["snakemake"] + argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            # Write output to temporary log file, to be attached later.
            log_file_path = os.path.join(config.log_path, "digestiflow-demux-snakemake.log.gz")
            with gzip.open(log_file_path, "wb") as log_file:
                mutex = Lock()
                logger_stderr = async_tee_pipe(proc, proc.stderr, log_file, sys.stderr, mutex)
                logger_stdout = async_tee_pipe(proc, proc.stdout, log_file, sys.stdout, mutex)
                logger_stderr.join()
                logger_stdout.join()
            # Copy out log file to log directory.
            os.makedirs(output_log_dir, exist_ok=True)
            shutil.copy(log_file_path, output_log_dir)
            failure = proc.returncode != 0
        except WorkflowError as e:
            logging.warning("Running demultiplexing failed: %s", e)
            failure = True

    if not failure and not config.api_read_only:
        message = send_flowcell_success_message(client, flowcell, output_dir, log_file_path)
        logging.info("Marking flowcell as complete...")
        try:
            client.flowcell_update(flowcell["sodar_uuid"], status_conversion="complete")
        except ApiException as e:
            logging.warning("Could not update conversion state to complete via API: %s", e)
        logging.info("Done running Snakemake.")
    elif flowcell and not config.api_read_only:
        message = send_flowcell_failure_message(client, flowcell, log_file_path)
        logging.info("Marking flowcell as failed...")
        try:
            client.flowcell_update(flowcell["sodar_uuid"], status_conversion="failed")
        except ApiException as e:
            logging.warning("Could not update conversion state to failed via API: %s", e)
    else:
        message = None
    return (not failure, message, flowcell, client)


def perform_demultiplexing(config, input_dir, output_dir):
    """Prepare and execute the demultiplexing with the Snakemake workflow."""
    logging.info("Starting to process input directory %s", input_dir)
    logging.info("Output will go to %s", output_dir)

    logging.debug("Creating output directory %s", output_dir)
    os.makedirs(output_dir, exist_ok=True)

    flowcell = create_sample_sheet(config, input_dir, output_dir)
    if not flowcell:
        return False, None, None, None

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
