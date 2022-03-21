"""Code for supporting the Snakemake file."""

import functools
import glob
import os

from .workflow import load_run_parameters

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Smallest RTA version to trigger call of bcl2fastq v2.
RTA_MIN_BCL2FASTQ2 = (1, 18, 54)


def listify(fn=None, wrapper=list):
    """
    A decorator which wraps a function's return value in ``list(...)``.

    Useful when an algorithm can be expressed more cleanly as a generator but
    the function should return an list.

    Example::

        >>> @listify
        ... def get_lengths(iterable):
        ...     for i in iterable:
        ...         yield len(i)
        >>> get_lengths(["spam", "eggs"])
        [4, 4]
        >>>
        >>> @listify(wrapper=tuple)
        ... def get_lengths_tuple(iterable):
        ...     for i in iterable:
        ...         yield len(i)
        >>> get_lengths_tuple(["foo", "bar"])
        (3, 3)
    """

    def listify_return(fn):
        @functools.wraps(fn)
        def listify_helper(*args, **kw):
            return wrapper(fn(*args, **kw))

        return listify_helper

    if fn is None:
        return listify_return
    return listify_return(fn)


def config_to_rta_version(config):
    """Given the flow cell configuration, return RTA version tuple."""
    input_dir = config["input_dir"]
    path_run_info = glob.glob(os.path.join(input_dir, "?un?arameters.xml"))[0]
    run_parameters = load_run_parameters(path_run_info)
    rta_version = run_parameters["rta_version"]
    return rta_version


def bcl2fastq_wrapper(config):
    """Return name of bcl2fastq wrapper to use."""
    rta_version = config_to_rta_version(config)
    if rta_version >= RTA_MIN_BCL2FASTQ2:
        return "bcl2fastq2"
    else:
        return "bcl2fastq"


def wrapper_path(path):
    """Generate path to wrapper"""
    return "file://" + os.path.abspath(os.path.join(os.path.dirname(__file__), "wrappers", path))


def undetermined_libraries(flowcell, rta_version):
    """Return library dicts for undetermined libraries in ``flowcell``."""
    lanes = set()
    for library in flowcell["libraries"]:
        lanes |= set(library["lanes"])
    result = []
    for lane in lanes:
        result.append(
            {
                "name": "lane{}".format(lane)
                if rta_version < RTA_MIN_BCL2FASTQ2
                else "Undetermined",
                "reference": library["reference"],
                "barcode": "Undetermined",
                "barcode2": "Undetermined",
                "lanes": [lane],
            }
        )
    return result


@listify
def lib_file_names(library, rta_version, n_template, n_index, lane=None, seq=None, name=None):
    """Return list with file names for the given library."""
    if rta_version < RTA_MIN_BCL2FASTQ2 and library.get("barcode2", "Undetermined") not in (
        "",
        "Undetermined",
    ):
        indices = ["".join((library["barcode"], "-", library["barcode2"]))]
    else:
        indices = [library["barcode"] or "NoIndex"]
    reads = ["R" + str(i + 1) for i in range(n_template)]
    reads += ["I" + str(i + 1) for i in range(n_index)]
    lanes = ["L{:03d}".format(lno) for lno in library["lanes"] if lane is None or lno == lane]
    if seq is None:
        seq = ""
    else:
        seq = seq + "_"
    if rta_version < RTA_MIN_BCL2FASTQ2:
        tpl = "{sample_name}_{index}_{lane}_{read}_001.fastq.gz"
    else:
        tpl = "{sample_name}_{seq}{lane}_{read}_001.fastq.gz"
    return list(
        sorted(
            [
                tpl.format(
                    sample_name=name or library["name"], index=index, lane=lane, read=read, seq=seq
                )
                for index in indices
                for read in reads
                for lane in lanes
            ]
        )
    )


def build_sample_map(flowcell):
    """Build sample map ``dict`` for the given flowcell."""
    result = {}
    rows = [(lane, lib["name"]) for lib in flowcell["libraries"] for lane in lib["lanes"]]
    i = 1
    for _, name in sorted(set(rows)):
        if name not in result:
            result[name] = "S{}".format(i)
            i += 1
    return result


@listify
def get_result_files_demux_with_project(config):
    """Return list with demultiplexing results."""

    def out_prefix(path):
        return os.path.join(config["output_dir"], path)

    flowcell = config["flowcell"]
    sample_map = build_sample_map(flowcell)
    bases_mask = flowcell["demux_reads"]
    n_template = bases_mask.count("T")
    n_index = (
        bases_mask.count("B")
        if (
            config["bcl2fastq2_params"]["create_fastq_for_index_reads"]
            and config["demux_tool"] == "bcl2fastq2"
        )
        else 0
    )
    expect_undetermined = True if "B" in bases_mask else False
    rta_version = config_to_rta_version(config)
    undetermined = undetermined_libraries(flowcell, rta_version) if expect_undetermined else []

    for lib in flowcell["libraries"] + undetermined:
        for lane in sorted(lib["lanes"]):
            if config["lanes"] and lane not in config["lanes"]:
                continue  # skip disabled lanes
            if lib["barcode"] == "Undetermined":
                sample_name = "Undetermined"
            else:
                sample_name = lib["name"]

            out_dir = (
                "{output_dir}/Undetermined/{id}".format(
                    id=flowcell["composite_id"],
                    output_dir=config["output_dir"],
                )
                if sample_name == "Undetermined"
                else "{output_dir}/{project}/{id}/{sample_name}".format(
                    id=flowcell["composite_id"],
                    project=lib.get("project_id", "Project"),
                    output_dir=config["output_dir"],
                    sample_name=sample_name,
                )
            )

            seq = sample_map.get(sample_name, "S0")
            name = "Undetermined" if lib["barcode"] == "Undetermined" else lib["name"]
            for fname in lib_file_names(lib, rta_version, n_template, n_index, lane, seq, name):
                yield out_prefix("{out_dir}/{fname}".format(out_dir=out_dir, fname=fname))


@listify
def get_result_files_demux(config):
    """Return list with demultiplexing results."""

    def out_prefix(path):
        return os.path.join(config["output_dir"], path)

    flowcell = config["flowcell"]
    sample_map = build_sample_map(flowcell)
    bases_mask = flowcell["demux_reads"]
    n_template = bases_mask.count("T")
    n_index = (
        bases_mask.count("B")
        if (
            config["bcl2fastq2_params"]["create_fastq_for_index_reads"]
            and config["demux_tool"] == "bcl2fastq2"
        )
        else 0
    )  # TODO check picard
    expect_undetermined = True if "B" in bases_mask else False
    rta_version = config_to_rta_version(config)
    undetermined = undetermined_libraries(flowcell, rta_version) if expect_undetermined else []

    for lib in flowcell["libraries"] + undetermined:
        for lane in sorted(lib["lanes"]):
            if config["lanes"] and lane not in config["lanes"]:
                continue  # skip disabled lanes
            if lib["barcode"] == "Undetermined":
                sample_name = "Undetermined"
            else:
                sample_name = lib["name"]
            out_dir = "{output_dir}/{sample_name}/{flowcell}/L{lane:03d}".format(
                flowcell=flowcell["vendor_id"],
                output_dir=config["output_dir"],
                sample_name=sample_name,
                lane=lane,
            )

            if rta_version < RTA_MIN_BCL2FASTQ2 or config["demux_tool"] == "picard":
                for fname in lib_file_names(lib, rta_version, n_template, n_index, lane):
                    yield out_prefix("{out_dir}/{fname}".format(out_dir=out_dir, fname=fname))
            else:
                seq = sample_map.get(sample_name, "S0")
                name = "Undetermined" if lib["barcode"] == "Undetermined" else lib["name"]
                for fname in lib_file_names(lib, rta_version, n_template, n_index, lane, seq, name):
                    yield out_prefix("{out_dir}/{fname}".format(out_dir=out_dir, fname=fname))


def get_result_files_fastqc(config):
    """Return list with FASTQC results files."""
    res_zip = []
    res_html = []
    paths = (
        get_result_files_demux(config)
        if not config["demux_tool"] == "bclconvert"
        else get_result_files_demux_with_project(config)
    )
    for path in paths:
        ext = ".fastq.gz"
        if path.endswith(ext):
            folder = os.path.dirname(path)
            base = os.path.basename(path)[: -len(ext)]
            res_zip.append(os.path.join(folder, "qc", "fastqc", base + "_fastqc.zip"))
            res_html.append(os.path.join(folder, "qc", "fastqc", base + "_fastqc.html"))
    return {"zip": res_zip, "html": res_html}


def get_tiles_arg(config):
    """Return tiles argument."""
    # Check whether all lanes are specified in the sheet
    lanes = set()
    for lib in config["flowcell"]["libraries"]:
        lanes |= set(lib["lanes"])
    all_lanes_in_sheet = len(lanes) == config["flowcell"]["num_lanes"]

    # Shortcut to whether any args are given for tiles
    args_given = bool(config["tiles"] or config["lanes"])

    # Handle case of no args given and fake --lane parameter if not all lanes mentioned in sheet.
    if not args_given and all_lanes_in_sheet:
        return ""  # no args
    elif not args_given and not all_lanes_in_sheet:
        config["lanes"] = list(lanes)
    # Handle case of selected lanes or selected tiles
    if config["tiles"]:
        return "--tiles {}".format(",".join(config["tiles"]))
    else:
        regexes = ["s_{}".format(lane) for lane in sorted(config["lanes"])]
        return "--tiles {}".format(",".join(regexes))


def get_tool_marker(config):
    """Return marker file for either bcl2fastq, bcl2fastq2, bclconvert or picard for snakemake"""
    marker_dict = {
        "bcl2fastq1": "bcl2fastq1.done",
        "bcl2fastq2": "bcl2fastq2.done",
        "bclconvert": "bclconvert.done",
        "picard": "picard.done",
    }
    return marker_dict[config["demux_tool"]]
